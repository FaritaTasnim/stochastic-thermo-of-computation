import tables as tb
import numpy as np
from scipy.stats import norm, binom, randint, uniform
from scipy.sparse.linalg import expm, expm_multiply, eigs
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix, csc_matrix
from scipy.sparse import identity, diags
from itertools import combinations, product, permutations
from copy import deepcopy
from tqdm import trange, tqdm
from functools import reduce
from time import time
from collections import namedtuple
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
import ternary
font = {'weight': 'bold',
        'size'   : 16}
plt.rc('font', **font)
# from numba import njit, jit
# %load_ext line_profiler
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

ValidNexts = namedtuple('ValidNexts', ('indices, states, propensities'))
Indices = namedtuple('Indices', ('ind', 'ind_io', 'ind_i', 'ind_o'))
Quant_Dots = namedtuple('Quant_Dots', ('q_dot', 'if_io_dot', 'if_oi_dot', 'mi_dot', 'mi_inst', 'zdot'))

rng = np.random.default_rng()

class Quants:
    def __init__(self, s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind):
        self.s = s
        self.q = q
        self.ep = ep
        self.if_io = if_io
        self.if_oi = if_oi
        self.mi_inst = mi_inst
        self.mi_tot = mi_tot
        self.z = z
        self.z_trans = z_trans
        self.nwind = nwind

class Intrxn_Struct():

    def joint_state_space (self, sp_spcs):
        '''sp_spcs: (dict) {spin: (nparray) that spin's state space for all spins}'''
        return np.array([np.array(s) for s in product(*sp_spcs.values())])

    def sub_state_index(self, sub, substate):
        return np.where(self.spin_state_spaces[sub]==substate)[0][0]

    def state_index(self, state):
        # print([np.where(self.X[:,i]==s)[0] for i,s in enumerate(state)])
        return reduce(np.intersect1d, [np.where(self.X[:,i]==s)[0] for i,s in enumerate(state)])[0]

    def state_indices(self, mstate, msubset, jstate_spc):
        if len(msubset) > 1:
            return reduce(np.intersect1d, [np.where(jstate_spc[:,sub]==mstate[i])[0] for i,sub in enumerate(msubset)])
        else:
            return reduce(np.intersect1d, [np.where(jstate_spc[:,sub]==mstate)[0] for i,sub in enumerate(msubset)])

    def get_marginal(self, subset, marg_spc, joint_spc, joint_pdist):
        return np.array([np.sum(np.take(joint_pdist, self.state_indices(ms, subset, joint_spc))) for ms in marg_spc])

    def set_mechs (self, iv):
        self.iv = iv # all (subsystem, mech_before, mech_after) triples where mechanism is represented by
                    # mech_before: a rule (float or function) controlling how rates are based on the before_state, and
                    # mech_after: a rule (float or function) controlling how rates are based on the after_state

    def where_spin_i_flips (self, spin_idx):
        # explanation: see https://colab.research.google.com/drive/1W-_NAn1JzyVRoc221qj4cRfS4WTUV-6p#scrollTo=nfXu6H21MXyZ
        s = self.spin_state_spaces
        s_lens = np.array([len(space) for space in s.values()]) # number of states each spin can take on

        num_diags = s_lens[spin_idx]-1 # number of positive-offset diagonals this spin will have
        offset_spacing = np.prod(s_lens[spin_idx+1:]) # spacing between each positive-offset diagonal
        offsets = np.array([offset_spacing*(i+1) for i in range(num_diags)]) # locations of positive offsets

        num_1_blocks = np.prod(s_lens[:spin_idx]) # number of 1 blocks filled with 1's
        num_blocks = num_1_blocks*2 - 1 # total number of blocks, because there will be num_blocks - 1 blocks filled with 0's

        def num_1s_and_0s_in_diag (this_diag): # gets the number of 1's and 0's in each 1 block and 0 block in a diagonal
            return (offsets[num_diags - 1 - this_diag], offsets[this_diag])

        diagonals = np.array([np.array([1 - b%2 for b in range(num_blocks)
                        for i in range(num_1s_and_0s_in_diag(d)[b%2])]) for d in range(num_diags)]) # construct diagonals

        # since the above was for positive-offset diagonals, and the negative-offset ones are completely symmetric,
        # the following doubles them up:
        if len(diagonals) == 1:
            all_diags = np.vstack((diagonals, diagonals))
        else:
            all_diags = np.hstack((diagonals, diagonals))

        all_offsets = np.hstack((offsets, -1*offsets))

        # use scipy.diags to construct the final M(i) matrix
        return diags(all_diags, all_offsets)


    def iv_diagonal (self, i, mech):
        if isinstance(mech, list):
            return [[diags(m_bef), diags(m_aft)] for m_bef, m_aft in mech]
        elif callable(mech):
            return diags([mech(i, state, self.intrxns['weighted']) for state in self.X])
        elif isinstance(mech, csc_matrix) or isinstance(mech, csr_matrix):
            return mech
        else:
            return diags(mech*np.ones(self.num_states))


    def set_all_RMs (self):
        Ks_iv = []
        for c, (i, v) in enumerate(self.iv):

            transitions_i = self.where_spin_i_flips(i)
            diags_iv = self.iv_diagonal(i, v)

            if isinstance(v, list):
                K_ivs = [d_bef.dot(transitions_i.dot(d_aft)) for d_bef, d_aft in diags_iv]
                # K_ivs = [d_aft.dot(transitions_i.dot(d_bef)) for d_bef, d_aft in diags_iv]
                K_iv = diags(np.zeros(self.num_states))
                for K_ in K_ivs:
                    K_iv += K_
                K_iv = K_iv + diags(np.array(-1*K_iv.sum(axis = 0))[0])
                Ks_iv.append(K_iv)

            elif isinstance(v, csc_matrix) or isinstance(v, csr_matrix):
                K_iv = v # blue must be normalized
                Ks_iv.append(K_iv)
            else:
                K_iv = transitions_i.dot(diags_iv)
                K_iv = K_iv + diags(np.array(-1*K_iv.sum(axis = 0))[0])
                Ks_iv.append(K_iv)

            if c == 0:
                K = K_iv
            else:
                K = K + K_iv

            c += 1

        self.Ks_iv = Ks_iv
        self.Ksys = K

    def normalized_j(self, b, a, K, p):
        j = K[a, b]*p[b] - K[b, a]*p[a]
        return j/p[a]

    def set_init_dist(self, id, p=None):
        # for pdist evolution you actually have to construct the pmf
        self.id = id

        if id != 'Cness' and id != 'test':
            init_dist = np.zeros(shape=(self.num_states, ))

            if id == 'uni-delta0':
                np.put(init_dist, np.where(self.X[:,self.o]==0)[0], 1.0)
                init_dist[self.invalid_state_inds] = 0.
            elif id == 'deltax0-delta0':
                np.put(init_dist, np.intersect1d(np.where(self.X[:,self.o]==0)[0], np.where(self.X[:,self.i]==self.Xi[1])[0]), 1.0)
            elif id == 'uni':
                init_dist = np.ones(shape=(self.num_states, ))
                init_dist[self.invalid_state_inds] = 0.
            elif id == 'delta0-delta0':
                np.put(init_dist, np.intersect1d(np.where(self.X[:,self.o]==0)[0], np.where(self.X[:,self.i]==0)[0]), 1.0)

            self.p_init = init_dist / np.sum(init_dist)

        else:
            self.p_init = p

    def sample_init_dist(self, init_dist = None):
        if init_dist is None:
            index = np.random.choice(self.num_states, p=self.p_init)
        else:
            index = np.random.choice(self.num_states, p=init_dist)
        state = self.X[index]
        return index, state

    def get_to_NESS (self, timestep, L1_thresh):
        # do pnew = p*e^(K*timestep) repeatedly until L1(pnew, p) < L1_thresh

        p = self.p_init
        t = 0
        L1 = 100
        self.L1_thresh = L1_thresh

        def step(p, t):
            p_bef = p
            p = expm_multiply((self.Ksys.tocsc()*timestep).tocsr(), p)
            p_aft = p
            L1 = sum([np.abs(x-y) for x,y in zip(p_bef, p_aft)])
            t += timestep
            return L1, p_aft, t

        while L1 > L1_thresh:
            L1, p, t = step(p, t)

        self.time_to_NESS = t
        self.p_ness = p
        self.ps_i_ness = {i: self.get_marginal([i], self.spin_state_spaces[i], self.X, self.p_ness) for i in range(self.num_spins)}
        self.p_io_ness = self.get_marginal(self.io, self.Xio, self.X, self.p_ness)

        # print(np.round(self.p_init,2))
        # print(np.round(self.p_ness,2))

        # self.norm_js = np.array([np.array([self.normalized_j(b, a, self.Ksys, p) for a in range(self.num_states)]) for b in range(self.num_states)])
        return t, self.p_ness, self.ps_i_ness, self.p_io_ness

    def get_state_tau_ago (self, traj, curr_t, tau):
        target = curr_t - tau
        return traj['states'][np.where(traj['times'] < target)[0][-1]]

    def row_append_ness(self, row, *args):
        if len(args) == 6:
            time, state, inds, iv, i, quants = args
            s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind = quants.s, quants.q, quants.ep, quants.if_io, quants.if_oi, quants.mi_inst, quants.mi_tot, quants.z, quants.z_trans, quants.nwind
            ind, ind_io, ind_i, ind_o = inds
        elif len(args) == 18:
            ep, if_io, if_oi, ind, ind_io, ind_i, ind_o, iv, mi_inst, mi_tot, nwind, q, s, state, time, i, z, z_trans = args
            # time, state, ind, ind_io, ind_i, ind_o, iv, i, s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind = args

        row['ind'] = ind
        row['ind_io'] = ind_io
        row['ind_i'] = ind_i
        row['ind_o'] = ind_o
        row['time'] = time
        row['state'] = state
        row['mech'] = iv
        row['who'] = i
        row['s'] = s
        row['q'] = q
        row['ep'] = ep
        row['if_io'] = if_io
        row['if_oi'] = if_oi
        row['mi_inst'] = mi_inst
        row['mi_tot'] = mi_tot
        row['z'] = z
        row['z_trans'] = z_trans
        row['nwind'] = nwind

        row.append()

    def row_append_rlx(self, row, *args):
        if len(args) == 10:
            time, state, inds, iv, i, quants, p, p_io, p_i, p_o = args
            s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind = quants.s, quants.q, quants.ep, quants.if_io, quants.if_oi, quants.mi_inst, quants.mi_tot, quants.z, quants.z_trans, quants.nwind
            ind, ind_io, ind_i, ind_o = inds
        elif len(args) == 22:
            ep, if_io, if_oi, ind, ind_io, ind_i, ind_o, iv, mi_inst, mi_tot, nwind, p, p_io, p_i, p_o, q, s, state, time, i, z, z_trans = args
            # time, state, ind, ind_io, ind_i, ind_o, iv, i, s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind, p, p_io, p_i, p_o = args

        row['ind'] = ind
        row['ind_io'] = ind_io
        row['ind_i'] = ind_i
        row['ind_o'] = ind_o
        row['time'] = time
        row['state'] = state
        row['mech'] = iv
        row['who'] = i
        row['s'] = s
        row['q'] = q
        row['ep'] = ep
        row['if_io'] = if_io
        row['if_oi'] = if_oi
        row['mi_inst'] = mi_inst
        row['mi_tot'] = mi_tot
        row['z'] = z
        row['z_trans'] = z_trans
        row['nwind'] = nwind
        row['p'] = p
        row['p_io'] = p_io
        row['p_i'] = p_i
        row['p_o'] = p_o

        row.append()

    def get_state_inds (self, state):
        joint = self.state_index(state)
        state_io = np.take(state, self.io)
        io = self.state_indices(state_io, [0,1], self.Xio)[0] # using 0,1 here because the state space Xio is for io already
        i = np.where(self.Xi==state[self.i])[0][0]
        o = np.where(self.Xo==state[self.o])[0][0]
        inds = Indices(joint, io, i, o)
        return inds

    def get_dots (self, p_start, curr_state, next_state, iv, flipper, sojourn, p=None, p_io=None, ps_i=None):

        if p is None:
            p = self.p_ness
            p_io = self.p_io_ness
            ps_i = self.ps_i_ness

        p_i, p_o = ps_i[self.i], ps_i[self.o]

        bef_inds = self.get_state_inds(curr_state)
        aft_inds = self.get_state_inds(next_state)
        bef, bef_io, bef_i, bef_o = bef_inds
        aft, aft_io, aft_i, aft_o = aft_inds

        s = np.log(p_start/p[aft])
        if iv == -1:
            q_dot = 0
        else:
            q_dot = np.log(self.Ks_iv[iv][aft,bef]/self.Ks_iv[iv][bef,aft])
            # if flipper == self.o:
            #     q_dot = (self.state_energies[bef] - self.state_energies[aft])/self.To[iv-1]
            # else:
            #     q_dot = 0
            # the above can be used for R1 and NR because the appropriate type of LDB holds in each
            # but leads to the same issue as the regular calculation of q_dot...
            # does that mean there's some error with the propensity scale when
            #  wondering if fundamental issue with q calculation for Gillespie when there's multiple heat baths per subsystem
        if_io_dot = 0
        if_oi_dot = 0

        zdot = np.log(p[bef] / p[aft]) + q_dot
        if flipper == self.o:
            if_io_dot = np.log((p_io[aft_io]*p_o[bef_o])/(p_io[bef_io]*p_o[aft_o]))
        if flipper == self.i:
            if_oi_dot = np.log((p_io[aft_io]*p_i[bef_i])/(p_i[aft_i]*p_io[bef_io]))
        mi_inst = np.log(p_io[bef_io]/(p_i[bef_i] * p_o[bef_o]))
        mi_dot = sojourn*mi_inst

        quant_dots = Quant_Dots(q_dot, if_io_dot, if_oi_dot, mi_dot, mi_inst, zdot)

        return bef_inds, aft_inds, s, quant_dots


    def get_options(self, ind):
        # print(np.round(self.Ksys.toarray(), 3))
        # all_propensities = np.concatenate([self.Ks_iv[i].getcol(ind).toarray().T[0] for i in range(len(self.iv))], axis=0)
        all_propensities = np.concatenate([K_iv.getcol(ind).toarray().T[0] for K_iv in self.Ks_iv], axis=0)
        valid_next_indices = np.where(all_propensities > 0)
        # valid_next_indices = np.mod(np.where(all_propensities > 0)[0], self.num_states)
        # valid_next_states = self.X[valid_next_indices]
        valid_next_states = self.X[np.mod(valid_next_indices[0], self.num_states)]
        valid_propensities = all_propensities[all_propensities > 0]
        sum_propensities = np.sum(valid_propensities)
        valid_nexts = ValidNexts(valid_next_indices, valid_next_states, valid_propensities)
        return all_propensities, valid_nexts, sum_propensities

    def get_sojourn(self, sum_propensities):
        r = rng.random()
        return -1 * np.log(r) / sum_propensities

    def get_transition(self, curr_state, valid_nexts, sum_propensities):
        valid_next_indices, valid_next_states, valid_propensities = valid_nexts
        propensity_scale = np.cumsum(valid_propensities)/sum_propensities
        # print(valid_next_indices)
        # print(propensity_scale)
        r = rng.random()
        choice = np.where(propensity_scale > r)[0].min()
        # print(choice)
        next_state = valid_next_states[choice]
        which_iv = valid_next_indices[0][choice] // self.num_states
        # print(which_iv)
        which_spin_flipped = np.where(curr_state != next_state)[0][0]
        return next_state, which_iv, which_spin_flipped

    def set_init_quants(self):
        s = 0
        q = 0
        ep = 0
        if_io = 0
        if_oi = 0
        mi_inst = 0
        mi_tot = 0
        z = np.zeros(self.num_spins)
        z_trans = np.zeros((self.num_states, self.num_states))
        nwind = 0

        init_quants = Quants(s, q, ep, if_io, if_oi, mi_inst, mi_tot, z, z_trans, nwind)
        return init_quants

    def set_traj_inits(self, p_init):
        init_index, init_state = self.sample_init_dist(p_init)
        p_init_state = p_init[init_index]

        init, init_io, init_i, init_o = self.get_state_inds(init_state)
        init_indices = Indices(init, init_io, init_i, init_o)

        t = 0
        curr_state = init_state
        next_state = None

        init_quants = self.set_init_quants()

        return p_init_state, init_indices, t, curr_state, next_state, init_quants

    def update_quants(self, quants, quant_dots, which_spin_flipped, a, b):
        quants.q += quant_dots.q_dot
        quants.ep = quants.s + quants.q
        quants.if_io += quant_dots.if_io_dot
        quants.if_oi += quant_dots.if_oi_dot
        quants.mi_inst = quant_dots.mi_inst
        quants.mi_tot += quant_dots.mi_dot
        quants.z[which_spin_flipped] += quant_dots.zdot
        quants.z_trans[a, b] += quant_dots.zdot
        if which_spin_flipped == self.i:
            quants.nwind += 1
        return quants


    ############################################################################


    def gentraj_NESS (self, row, table):

        p_init_state, init_indices, t, curr_state, next_state, quants = self.set_traj_inits(self.p_ness)
        self.row_append_ness(row, t, curr_state, init_indices, -2, -2, quants)

        while t < self.ness_runtime:
            if t > 0:
                curr_state = next_state

            curr_ind = self.state_index(curr_state)
            all_propensities, valid_nexts, sum_propensities = self.get_options(curr_ind)

            sojourn = self.get_sojourn(sum_propensities)
            t += sojourn

            if t < self.ness_runtime:
                next_state, which_iv, which_spin_flipped = self.get_transition(curr_state, valid_nexts, sum_propensities)
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, next_state, which_iv, which_spin_flipped, sojourn)
                quants = self.update_quants(quants, quant_dots, which_spin_flipped, aft.ind, bef.ind)
                self.row_append_ness(row, t, next_state, aft, which_iv, which_spin_flipped, quants)

            else:
                soj = self.ness_runtime - t + sojourn
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, curr_state, -1, -1, soj)
                quants.mi_inst = quant_dots.mi_inst
                quants.mi_tot += quant_dots.mi_dot
                self.row_append_ness(row, self.ness_runtime, curr_state, bef, -2, -2, quants)

    def genwind_NESS (self, row, table):

        p_init_state, init_indices, t, curr_state, next_state, quants = self.set_traj_inits(self.p_ness)
        self.row_append_ness(row, t, curr_state, init_indices, -2, -2, quants)

        while t < 1000000:
            if t > 0:
                curr_state = next_state

            curr_ind = self.state_index(curr_state)
            all_propensities, valid_nexts, sum_propensities = self.get_options(curr_ind)

            sojourn = self.get_sojourn(sum_propensities)
            t += sojourn

            next_state, which_iv, which_spin_flipped = self.get_transition(curr_state, valid_nexts, sum_propensities)

            if which_spin_flipped == self.i:
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, curr_state, -1, -1, sojourn)
                quants.mi_inst = quant_dots.mi_inst
                quants.mi_tot += quant_dots.mi_dot
                self.row_append_ness(row, self.ness_runtime, curr_state, bef, -2, -2, quants)
                break
            else:
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, next_state, which_iv, which_spin_flipped, sojourn)
                quants = self.update_quants(quants, quant_dots, which_spin_flipped, aft.ind, bef.ind)
                self.row_append_ness(row, t, next_state, aft, which_iv, which_spin_flipped, quants)

    def row_append_pdist_relax(self, row, time, pdist, pdist_io, pdist_i, pdist_o):
        row['time'] = time
        row['p'] = pdist
        row['p_io'] = pdist_io
        row['p_i'] = pdist_i
        row['p_o'] = pdist_o
        row.append()

    def generate_pdists_relaxation(self, row, table, num_steps):
        p = self.p_init
        p_io = self.get_marginal(self.io, self.Xio, self.X, p)
        p_i = self.get_marginal([self.i], self.Xi , self.X, p)
        p_o = self.get_marginal([self.o], self.Xo , self.X, p)
        t = 0
        self.row_append_pdist_relax(row, t, p, p_io, p_i, p_o)
        timestep = self.time_to_NESS / num_steps

        def step(p, t):
            p_bef = p
            p = expm_multiply((self.Ksys.tocsc()*timestep).tocsr(), p)
            p_aft = p
            p_aftio = self.get_marginal(self.io, self.Xio, self.X, p)
            p_afti = self.get_marginal([self.i], self.Xi , self.X, p)
            p_afto = self.get_marginal([self.o], self.Xo , self.X, p)
            L1 = sum([np.abs(x-y) for x,y in zip(p_bef, p_aft)])
            t += timestep
            return L1, p_aft, p_aftio, p_afti, p_afto, t

        while t <= self.time_to_NESS:
            L1, p, p_io, p_i, p_o, t = step(p, t)
            self.row_append_pdist_relax(row, t, p, p_io, p_i, p_o)

    def get_all_ps(self, p):
        p_io = self.get_marginal(self.io, self.Xio, self.X, p)
        ps_i = {i: self.get_marginal([i], self.spin_state_spaces[i], self.X, p) for i in range(self.num_spins)}
        p_i = ps_i[self.i]
        p_o = ps_i[self.o]

        return p_io, ps_i, p_i, p_o

    def get_traj_rlx_inits(self):
        init_index, init_state = self.sample_init_dist()
        p = self.p_init
        p_init_state = p[init_index]
        p_io, ps_i, p_i, p_o = self.get_all_ps(self.p_init)
        init_inds = self.get_state_inds(init_state)

        t = 0
        quants = self.set_init_quants()

        return t, init_state, init_inds, quants, p_init_state, p, p_io, ps_i, p_i, p_o

    def gentraj_relax (self, row, table, saved_ts, saved_ps):

        t, curr_state, init_inds, quants, p_init_state, p, p_io, ps_i, p_i, p_o = self.get_traj_rlx_inits()
        self.row_append_rlx(row, t, curr_state, init_inds, -2, -2, quants, p, p_io, p_i, p_o)

        while t < self.time_to_NESS:
            if t > 0:
                curr_state = next_state

            curr_ind = self.state_index(curr_state)
            all_propensities, valid_nexts, sum_propensities = self.get_options(curr_ind)
            sojourn = self.get_sojourn(sum_propensities)
            t += sojourn

            nearest_ind = np.where(saved_ts <= t)[0][-1]
            extra_bitof_t = t - saved_ts[nearest_ind]
            p = saved_ps[nearest_ind]
            p = expm_multiply((self.Ksys.tocsc()*extra_bitof_t).tocsr(), p)
            # p = expm_multiply((self.Ksys.tocsc()*sojourn).tocsr(), p)
            p_io, ps_i, p_i, p_o = self.get_all_ps(p)

            if t < self.time_to_NESS:
                next_state, which_iv, which_spin_flipped = self.get_transition(curr_state, valid_nexts, sum_propensities)
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, next_state, which_iv, which_spin_flipped, sojourn, p, p_io, ps_i)
                quants = self.update_quants(quants, quant_dots, which_spin_flipped, aft.ind, bef.ind)
                self.row_append_rlx(row, t, next_state, aft, which_iv, which_spin_flipped, quants, p, p_io, p_i, p_o)

            else:
                soj = self.time_to_NESS - t + sojourn
                bef, aft, s, quant_dots = self.get_dots(p_init_state, curr_state, curr_state, -1, -1, soj, p, p_io, ps_i)
                quants.mi_inst = quant_dots.mi_inst
                quants.mi_tot += quant_dots.mi_dot
                self.row_append_rlx(row, self.time_to_NESS, curr_state, bef, -2, -2, quants, p, p_io, p_i, p_o)



class Lattice(Intrxn_Struct):

    def __init__(self, lattice_shape, spin_states, toroidal, input_spin, output_spin, energy_func, ness_runtime):
        # choose dimension of lattice, and side length of lattice, options are the set of positive integers, Z+
        # note this means that total number of spins is the product of all the dimensions of lattice_shape
        self.i = input_spin
        self.o = output_spin
        self.lattice_shape = lattice_shape
        self.num_spins = np.prod(self.lattice_shape)

        # generate interaction network
        self.system = [i for i in range(self.num_spins)]
        self.intrxns = self.generate_lattice(lattice_shape, toroidal)
        # get all valid units in the system, except the system itself
        # self.units = self.get_units(self.intrxns['unweighted'])
        # generate the state spaces for each spin, the whole system, and each unit
        self.num_spin_states = len(spin_states)
        self.spin_state_spaces = {i: np.array(spin_states[i]) for i in spin_states}
        self.X = self.joint_state_space(self.spin_state_spaces)
        self.num_states = len(self.X)

        self.Xi = self.spin_state_spaces[self.i]
        self.Xo = self.spin_state_spaces[self.o]
        self.num_istates = len(self.Xi)
        self.num_ostates = len(self.Xo)

        self.io = [self.i, self.o]
        self.Xio = self.joint_state_space({i: self.spin_state_spaces[i] for i in self.io})
        self.num_iostates = len(self.Xio)
        self.ness_runtime = ness_runtime

        if energy_func == 'abs':
            self.H = lambda state: np.abs(state[self.i] - state[self.o])
            self.Hname = 'abs'
        elif energy_func == 'harmonic':
            self.H = lambda state: (state[self.i] - state[self.o])**2
            self.Hname = 'harmonic'
        elif energy_func == 'delta':
            self.H = lambda state: -1 if state[self.i] == state[self.o] else 0
            self.Hname = 'delta'
        elif energy_func == 'nobias':
            self.H = lambda state: 0
            self.Hname = 'nobias'
        else:
            self.H = energy_func[1] # input a state, outputs an energy value
            self.Hname = energy_func[0]

        self.state_energies = np.apply_along_axis(self.H, 1, self.X)

    def generate_lattice (self, lat_shape, toroidal):
        ndim = len(lat_shape)

        # create a numpy ndarray to represent the regular lattice
        # allows for efficient operations over the lattice space
        lattice_array = np.arange(self.num_spins).reshape(lat_shape)

        # function to get the index (location) of a desired spin in the lattice
        loc = lambda s : list(zip(*np.where(lattice_array == s)))[0]

        all_spin_locs = [loc(i) for i in self.system]

        offsets = []
        for dim in range(ndim):
            offsetplus = np.zeros(ndim)
            offsetplus[dim] = 1
            offsets.append(offsetplus)
            offsetminus = np.zeros(ndim)
            offsetminus[dim] = -1
            offsets.append(offsetminus)
        offsets = np.array(offsets, dtype=int)

        def get_leaders(spin, toroidal):

            this_spin_loc = loc(spin)

            def correct_edges (row):
                out = deepcopy(row)

                for idx in range(len(self.lattice_shape)):
                    if row[idx] == self.lattice_shape[idx] and this_spin_loc[idx] == self.lattice_shape[idx] - 1:
                        out[idx] = 0

                return out

            if spin == self.i:
                leaders = []

            else:
                leaders = this_spin_loc + offsets
                if toroidal:
                    leaders = np.apply_along_axis(correct_edges, 1, leaders)

                else:
                    valid = np.all((leaders < np.array(lat_shape)) & (leaders >= 0), axis=1)
                    leaders = leaders[valid]

            return set([lattice_array.item(tuple(n)) for n in leaders])

        all_leaders = {i: get_leaders(i, toroidal) for i in self.system}

        return {'weighted': {i: {j: 1 for j in all_leaders[i]} for i in self.system},
                'unweighted': all_leaders}

class Chain(Lattice):

    def generate_lattice(self, num_spins):
        all_leaders = {i: [i-1, i] if i != 0 else [i] for i in self.system}
        return {'weighted': {i: {j: 1 for j in all_leaders[i]} for i in self.system},
                'unweighted': all_leaders}

class Comm_Channel(Lattice):

    def set_attrs(self, attrs):
        if 'sr' in attrs.keys():
            self.sr = attrs['sr']
        else:
            self.sr = 0
        if 'To' in attrs.keys():
            self.To = attrs['To']

    def set_group (self, f, root, typ):
        group = f.create_group(root, typ)
        group._v_attrs.space = self.X
        group._v_attrs.space_io = self.Xio
        group._v_attrs.space_i = self.Xi
        group._v_attrs.space_o = self.Xo
        group._v_attrs.p_io = self.p_io_ness
        group._v_attrs.p_i = self.ps_i_ness[self.i]
        group._v_attrs.p_o = self.ps_i_ness[self.o]
        group._v_attrs.p_init = self.p_init
        group._v_attrs.K = self.Ksys
        group._v_attrs.sr = self.sr
        group._v_attrs.To = self.To
        group._v_attrs.H = self.Hname
        group._v_attrs.state_energies = self.state_energies
        group._v_attrs.Cness = self.Cness
        return group

    def set_cols(self):
        self.trans_shape = (self.num_states, self.num_states)

        class cols_ness(tb.IsDescription):
            ep = tb.Float64Col()
            if_io = tb.Float64Col()
            if_oi = tb.Float64Col()
            ind = tb.Int16Col()
            ind_io = tb.Int16Col()
            ind_i = tb.Int16Col()
            ind_o = tb.Int16Col()
            mech = tb.Int8Col()
            mi_inst = tb.Float64Col()
            mi_tot = tb.Float64Col()
            nwind = ind = tb.Int16Col()
            q = tb.Float64Col()
            s = tb.Float64Col()
            state = tb.Int8Col(shape = (1, self.num_spins))
            time = tb.Float64Col()
            who = tb.Int8Col()
            z = tb.Float64Col(shape = (1, self.num_spins))
            z_trans = tb.Float64Col(shape = self.trans_shape)

        class cols_rlx(tb.IsDescription):
            ep = tb.Float64Col()
            if_io = tb.Float64Col()
            if_oi = tb.Float64Col()
            ind = tb.Int16Col()
            ind_io = tb.Int16Col()
            ind_i = tb.Int16Col()
            ind_o = tb.Int16Col()
            mech = tb.Int8Col()
            mi_inst = tb.Float64Col()
            mi_tot = tb.Float64Col()
            nwind = ind = tb.Int16Col()
            p = tb.Float64Col(shape = (1, self.num_states))
            p_io = tb.Float64Col(shape = (1, self.num_iostates))
            p_i = tb.Float64Col(shape = (1, self.num_istates))
            p_o = tb.Float64Col(shape = (1, self.num_ostates))
            q = tb.Float64Col()
            s = tb.Float64Col()
            state = tb.Int8Col(shape = (1, self.num_spins))
            time = tb.Float64Col()
            who = tb.Int8Col()
            z = tb.Float64Col(shape = (1, self.num_spins))
            z_trans = tb.Float64Col(shape = self.trans_shape)

        class cols_pdist(tb.IsDescription):
            p = tb.Float64Col(shape = (1, self.num_states))
            p_io = tb.Float64Col(shape = (1, self.num_iostates))
            p_i = tb.Float64Col(shape = (1, self.num_istates))
            p_o = tb.Float64Col(shape = (1, self.num_ostates))
            time = tb.Float64Col()

        class cols_rlx_qs(tb.IsDescription):
            ep = tb.Float64Col()
            if_io = tb.Float64Col()
            if_oi = tb.Float64Col()
            mi_inst = tb.Float64Col()
            mi_tot = tb.Float64Col()
            nwind = tb.Int16Col()
            q = tb.Float64Col()
            s = tb.Float64Col()
            time = tb.Float64Col()

        self.cols_ness = cols_ness
        self.cols_rlx = cols_rlx
        self.cols_pdist = cols_pdist
        self.cols_rlx_qs = cols_rlx_qs

    def write_pdists(self, fname, cols, cond, num_steps):
        tb.file._open_files.close_all()
        f = tb.open_file(fname, mode="w")
        root = f.root
        grp = self.set_group(f, root, 'pdists_relax')
        table = f.create_table(grp, 'pdists_relax', cols, 'pdists_relax')
        r = table.row
        self.generate_pdists_relaxation(r, table, num_steps)
        table.flush()
        f.close()

    def write_traj(self, fname, WorT, NorR, cols, cond, num_runs, saved_ts=None, saved_ps=None):
        '''
            model: (obj) instance of communication channel to be gathering data for
            filename: (str) name of file to create and save data to
            WorT: (str) either 'w' or 't' to indicate that windows or trajectories are to be collected
            NorR: (str) eithr 'n' or 'r' to indicate either NESS or relaxation process to be run
            cols: (obj) class of type tb.IsDescription that defines the columns for the trajectory data
            cond: (str) condition number from model_conditions dictionary
            group_fcn: (fcn) for setting the group containing all the trajectories
            group_attrs: (list) of attributes of the group
            get_traj_fcn: (fcn) for getting a single trajectory
            num_runs: (int) number of trajectories to generate
        '''
        def run(runfile, root, group, run_fcn, *args):
            for i in trange(num_runs):
                this_run = 'run_' + str(i).zfill(7)
                table = runfile.create_table(group, this_run, cols, this_run)
                r = table.row
                run_fcn(r, table, *args)
                table.flush()


        f = tb.open_file(fname, mode="w")
        root = f.root
        if WorT == 't':
            grp = self.set_group(f, root, 'runs')
            if NorR == 'n':
                run(f, root, grp, self.gentraj_NESS)
            elif NorR == 'r':
                run(f, root, grp, self.gentraj_relax, saved_ts, saved_ps)
        elif WorT == 'w':
            grp = self.set_group(f, root, 'runs')
            if NorR == 'n':
                run(f, root, grp, self.genwind_NESS)
            elif NorR == 'r':
                pass
        f.close()

    def write_NESS_qs (self, c, TorW):

        if TorW == 't':
            rf = tb.open_file(c + '_traj_NESS.h5', mode = "r")
            wf = tb.open_file(c + '_traj_quants_NESS.h5', mode = "w")
        elif TorW == 'w':
            rf = tb.open_file(c + '_wind_NESS.h5', mode = "r")
            wf = tb.open_file(c + '_wind_quants_NESS.h5', mode = "w")

        rfr = rf.root
        wfr = wf.root

        for attr in rfr.runs._v_attrs._f_list():
            wfr._v_attrs[attr] = rfr.runs._v_attrs[attr]
        wfg = wf.create_group(wfr, 'quants')
        wt = wf.create_table(wfg, 'quants', self.cols_ness, 'quants')
        for run in tqdm(rfr.runs):
            r = wt.row
            self.row_append_ness(r, *run[-1])
        wt.flush()

        rf.close()
        wf.close()


    def write_rlx_qs (self, c):

        pf = tb.open_file(c + '_pdists_relax.h5', mode='r')
        pr = pf.root

        tf = tb.open_file(c + '_traj_relax.h5', mode='r')
        tr = tf.root

        wf = tb.open_file(c + '_quants_relax.h5', mode='w')
        wfr = wf.root

        def get_quants(traj, quant, indices):
            vals = traj.col(quant)
            taken_vals = np.take(vals, indices)
            return taken_vals

        def relax_row_append(r, ep, if_io, if_oi, mi_inst, mi_tot, nwind, q, s, t):
            r['ep'] =  ep
            r['if_io'] = if_io
            r['if_oi'] = if_oi
            r['mi_inst'] = mi_inst
            r['mi_tot'] = mi_tot
            r['nwind'] = nwind
            r['q'] = q
            r['s'] = s
            r['time'] = t

            r.append()

        wfg = wf.create_group(wfr, 'quants')
        for attr in tr.runs._v_attrs._f_list():
            wfr._v_attrs[attr] = tr.runs._v_attrs[attr]

        timesteps = pr.pdists_relax.pdists_relax.col('time')
        pdists_at_times = pr.pdists_relax.pdists_relax.col('p')[:,0]
        piodists_at_times = pr.pdists_relax.pdists_relax.col('p_io')[:,0]
        pidists_at_times = pr.pdists_relax.pdists_relax.col('p_i')[:,0]
        podists_at_times = pr.pdists_relax.pdists_relax.col('p_o')[:,0]

        i = 0

        def to_full(a):
            output = np.full((len(a), max(map(len, a))), np.nan)
            for i, row in enumerate(a):
                output[i, :len(row)] = row
            return output

        for traj in tqdm(tr.runs):

            trans_times = traj.col('time')
            states = traj.col('state')[:,0]
            pdists_at_trans = traj.col('p')[:,0]
            piodists_at_trans = traj.col('p_io')[:,0]
            pidists_at_trans = traj.col('p_i')[:,0]
            podists_at_trans = traj.col('p_o')[:,0]
            inds = traj.col('ind')
            inds_io = traj.col('ind_io')
            inds_i = traj.col('ind_i')
            inds_o = traj.col('ind_o')

            num_steps = len(timesteps)

            # v below is intended to be a vector of [time1, time2]
            get_trans = lambda v: np.intersect1d(np.where(v[0] < trans_times)[0], np.where(trans_times < v[1])[0])
            timesteps_shift = np.concatenate([timesteps[1:], [timesteps[-1]]])
            timestep_pairs = np.concatenate([[timesteps], [timesteps_shift]], axis = 0).T
            all_trans_inds = to_full([get_trans(pair) for pair in timestep_pairs])

            nearest_trans = lambda x: np.where(trans_times <= x)[0][-1]
            vnt = np.vectorize(nearest_trans)
            nearest_trans_indices = vnt(timesteps_shift)

            num_intervals = all_trans_inds.shape[1] + 1

            wheretp = np.argwhere(~np.isnan(all_trans_inds))
            whattp = all_trans_inds[~np.isnan(all_trans_inds)]

            inds_io_intvl = [np.take(inds_io, nearest_trans_indices).reshape((num_steps, 1))]
            inds_i_intvl = [np.take(inds_i, nearest_trans_indices).reshape((num_steps, 1))]
            inds_o_intvl = [np.take(inds_o, nearest_trans_indices).reshape((num_steps, 1))]
            pdists_io_intvl = [piodists_at_times]
            pdists_i_intvl = [pidists_at_times]
            pdists_o_intvl = [podists_at_times]

            interval_caps = np.repeat([timesteps], num_intervals + 1, axis=0)
            interval_caps[-1] = timesteps_shift

            for intvl in range(all_trans_inds.shape[1]):

                nti = deepcopy(nearest_trans_indices)
                locs = wheretp[np.where(wheretp[:,1]==intvl)][:,0]
                vals = whattp[np.where(wheretp[:,1]==intvl)].astype(int)
                np.put(nti, locs, vals)

                t_intvl = intvl + 1
                interval_caps[t_intvl] = interval_caps[t_intvl-1]
                interval_caps[t_intvl, locs] = trans_times[vals]

                inds_io_intvl.append(np.take(inds_io, nti).reshape((num_steps, 1)))
                inds_i_intvl.append(np.take(inds_i, nti).reshape((num_steps, 1)))
                inds_o_intvl.append(np.take(inds_o, nti).reshape((num_steps, 1)))

                piodists = deepcopy(pdists_io_intvl[-1])
                pidists = deepcopy(pdists_i_intvl[-1])
                podists = deepcopy(pdists_o_intvl[-1])

                pio_replacers = np.take(piodists_at_trans, vals, axis=0)
                pi_replacers = np.take(pidists_at_trans, vals, axis=0)
                po_replacers = np.take(podists_at_trans, vals, axis=0)
                piodists[locs,:] = pio_replacers
                pidists[locs,:] = pi_replacers
                podists[locs,:] = po_replacers

                pdists_io_intvl.append(piodists)
                pdists_i_intvl.append(pidists)
                pdists_o_intvl.append(podists)

            time_intervals = (interval_caps[1:] - interval_caps[:-1]).T

            mi_vals_intvl = []
            for indices_io, indices_i, indices_o, pdists_io, pdists_i, pdists_o in zip(inds_io_intvl, inds_i_intvl, inds_o_intvl, pdists_io_intvl, pdists_i_intvl, pdists_o_intvl):
                pio_vals = np.take_along_axis(pdists_io, indices_io, axis=1)
                pi_vals = np.take_along_axis(pdists_i, indices_i, axis=1)
                po_vals = np.take_along_axis(pdists_o, indices_o, axis=1)
                out = np.zeros(pio_vals.shape)
                out2 = out = np.zeros(pio_vals.shape)
                pipo_vals = pi_vals * po_vals
                insidelog = np.divide(pio_vals, pipo_vals, out = out, where = pipo_vals!=0)
                mi_vals = np.log(insidelog, out = out2, where = insidelog!=0)
                # mi_vals = np.log(pio_vals / (pi_vals * po_vals))
                mi_vals_intvl.append(mi_vals)

            mis_intervals = np.concatenate(mi_vals_intvl, axis=1)
            mi_insts = mis_intervals[:,0]
            out = np.sum(mis_intervals*time_intervals, axis=1)
            mi_tots = np.cumsum(out)

            joint_inds = np.take(inds, nearest_trans_indices).reshape((num_steps, 1))
            p_init_state = pdists_at_times[0][inds[0]]
            ss = -1*np.log(np.take_along_axis(pdists_at_times, joint_inds, axis=1) / p_init_state).T

            if_ios = get_quants(traj, 'if_io', nearest_trans_indices)
            if_ois = get_quants(traj, 'if_oi', nearest_trans_indices)
            qs = get_quants(traj, 'q', nearest_trans_indices)
            nwinds = get_quants(traj, 'nwind', nearest_trans_indices)

            if_ios = if_ios.reshape((1, len(if_ios)))
            if_ois = if_ois.reshape((1, len(if_ois)))
            mi_tots = mi_tots.reshape(1, len(mi_tots))
            mi_insts = mi_insts.reshape(1, len(mi_insts))
            qs = qs.reshape((1, len(qs)))
            eps = qs + ss
            nwinds = nwinds.reshape((1, len(nwinds)))
            times = timesteps.reshape((1, len(timesteps)))

            quants = np.concatenate([eps, if_ios, if_ois, mi_insts, mi_tots, nwinds, qs, ss, times], axis=0).T

            this_run = 'run_' + str(i).zfill(7)
            wftable = wf.create_table(wfg, this_run, self.cols_rlx_qs, this_run)
            for item in quants:
                wfrow = wftable.row
                relax_row_append(wfrow, *item)
            wftable.flush()

            i += 1

        pf.close()
        tf.close()
        wf.close()

    def write_data(self, c, num_runs = 0, num_steps = 0, NESS = False, relax = False, wind = False, pdists = False, saved_ts = None, saved_ps = None):

        self.set_cols()
        if NESS:
            self.write_traj(c + '_traj_NESS.h5', 't', 'n', self.cols_ness, c, num_runs)
            self.write_NESS_qs(c, 't')
        if relax:
            self.write_traj(c + '_traj_relax.h5', 't', 'r', self.cols_rlx, c, num_runs, saved_ts, saved_ps)
            self.write_rlx_qs(c)
        if wind:
            self.write_traj(c + '_wind_NESS.h5', 'w', 'n', self.cols_ness, c, num_runs)
            self.write_NESS_qs(c, 'w')
        if pdists:
            self.write_pdists(c + '_pdists_relax.h5', self.cols_pdist, c, num_steps)

    def sanity_check(self, c):

        pf = tb.open_file(c + '_pdists_relax.h5', mode='r')
        pr = pf.root

        qf = tb.open_file(c + '_quants_relax.h5', mode='r')
        qr = qf.root

        def ensemble_quants_dots(p, p_io, p_i, p_o):

            def Ktrans_iv(bef, aft, iv):
                return self.Ks_iv[iv][aft, bef]

            def Jtrans_iv(bef, aft, iv):
                return Ktrans_iv(bef, aft, iv)*p[bef] - Ktrans_iv(aft, bef, iv)*p[aft]

            def Ktrans(bef, aft):
                return self.Ksys[aft, bef]

            def Jtrans(bef,aft):
                return Ktrans(bef, aft)*p[bef] - Ktrans(aft, bef)*p[aft]

            # the below methods assume the system consists only of two subsystems
            # particularly the calculation of ep and q with the for loop over the product space over states
            # will have to modify for general system
            if_io_dot = 0
            if_oi_dot = 0
            q_dot = 0
            ep_dot = 0
            mi_dot = 0

            for xi, xo, xpo in product(self.Xi, self.Xo, self.Xo):
                if xo != xpo:
                    bef = self.state_index(np.array([xi, xpo]))
                    bef_i = self.sub_state_index(self.i, xi)
                    bef_o = self.sub_state_index(self.o, xpo)
                    aft = self.state_index(np.array([xi, xo]))
                    aft_i = bef_i
                    aft_o = self.sub_state_index(self.o, xo)

                    if Ktrans(bef,aft) !=0 and p[bef] != 0 and p[aft] != 0:
                        pf = Jtrans(bef,aft)
                        if_io_dot += 0.5*pf*(np.log(p_io[aft]/p_o[aft_o]) - np.log(p_io[bef]/p_o[bef_o]))
                        pchange = np.log(p[bef]/p[aft])
                        for iv in range(len(self.iv)):
                            if Ktrans_iv(bef,aft,iv) != 0:
                                Kchange = np.log(Ktrans_iv(bef,aft,iv)/Ktrans_iv(aft,bef,iv))
                                pf = Jtrans_iv(bef,aft,iv)
                                q_dot += 0.5*pf*Kchange
                                ep_dot += 0.5*pf*(Kchange + pchange)

            for xpi, xi, xo in product(self.Xi, self.Xi, self.Xo):
                if xi != xpi:
                    bef = self.state_index(np.array([xpi, xo]))
                    bef_i = self.sub_state_index(self.i, xpi)
                    bef_o = self.sub_state_index(self.o, xo)
                    aft = self.state_index(np.array([xi, xo]))
                    aft_i = self.sub_state_index(self.i, xi)
                    aft_o = bef_o

                    if Ktrans(bef,aft) !=0 and p[bef] != 0 and p[aft] != 0:
                        pf = Jtrans(bef,aft)
                        if_oi_dot += 0.5*pf*(np.log(p_io[aft]/p_i[aft_i]) - np.log(p_io[bef]/p_i[bef_i]))
                        pchange = np.log(p[bef]/p[aft])
                        for iv in range(len(self.iv)):
                            if Ktrans_iv(bef,aft,iv) != 0:
                                Kchange = np.log(Ktrans_iv(bef,aft,iv)/Ktrans_iv(aft,bef,iv))
                                pf = Jtrans_iv(bef,aft,iv)
                                q_dot += 0.5*pf*Kchange
                                ep_dot += 0.5*pf*(Kchange + pchange)

            for (xi, xo) in self.X:
                ind = self.state_index(np.array([xi, xo]))
                ind_i = self.sub_state_index(self.i, xi)
                ind_o = self.sub_state_index(self.o, xo)

                if p_io[ind] !=0 and p_i[ind_i] != 0 and p_o[ind_o] != 0:
                    mi_dot += p_io[ind] * np.log(p_io[ind]/ (p_i[ind_i] * p_o[ind_o]))

            return np.array([ep_dot, if_io_dot, if_oi_dot, 0, mi_dot, 0, q_dot, 0])

        timesteps = pr.pdists_relax.pdists_relax.col('time')
        pdists_at_times = pr.pdists_relax.pdists_relax.col('p')[:,0]
        pidists_at_times = pr.pdists_relax.pdists_relax.col('p_i')[:,0]
        podists_at_times = pr.pdists_relax.pdists_relax.col('p_o')[:,0]

        step = timesteps[1]
        # print(timesteps)
        ensemble_quants_at_times = [np.array([0,0,0,0,0,0,0,0])]

        for i, (time, p, p_i, p_o) in tqdm(enumerate(zip(timesteps, pdists_at_times, pidists_at_times, podists_at_times))):
            ensemble_quants_at_times.append(ensemble_quants_at_times[-1] + step*ensemble_quants_dots(p, p, p_i, p_o))

        quants = ['ep', 'if_io', 'if_oi', 'mi_inst', 'mi_tot', 'nwind', 'q', 's']
        # each q is a table corresponding to one trajectory
        all = np.concatenate([[np.concatenate([[q.col(quant)] for quant in quants], axis = 0).T] for q in qr.quants], axis = 0)
        # shape (num_traj, num_timesteps, num_quants)

        ensemble_quants_at_times = np.array(ensemble_quants_at_times)[1:]

        s = np.array([sum(-1*pio*np.log(pio) if pio > 0 else 0 for pio in p_io) for p_io in pdists_at_times])
        sdot = s-s[0]
        ensemble_quants_at_times[:,7] = sdot

        ensemble_quants_at_times[:,3] = np.concatenate([[0],np.diff(ensemble_quants_at_times[:,4])/step])

        avgd_traj_quants_at_times = np.mean(all, axis=0) #average at each timestep (average over trajectories)

        ########################################################################

        for i in range(ensemble_quants_at_times.shape[1]):
            print(quants[i])
            print(np.round(ensemble_quants_at_times[:,i], 3))
            print(np.round(avgd_traj_quants_at_times[:,i], 3))
            print()

        print()

    def append_output_mechs(self, mechs):

        for To in self.To:

            aft = np.exp(0.5*self.state_energies/To)
            bef = np.exp(-0.5*self.state_energies/To)

            mechs.append( (self.o, [(bef, aft)]) )

        return mechs

    def expand_marginal(self, ss_i, p_i, which_sub):
        out = np.zeros(self.num_states)
        for ind, ival in enumerate(ss_i):
            x = np.where(self.X[:,which_sub] == ival)
            prob = p_i[ind]
            np.put(out, x, prob)
        return out

    def I (self, p_io=None, p_i=None, p_o=None):
        if p_io is None:
            p_io = self.p_io_ness
            p_i, p_o = self.ps_i_ness[self.i], self.ps_i_ness[self.o]
        mi = 0
        for idx, i in enumerate(self.Xi):
            this_pi = p_i[idx]
            if this_pi  > 0:
                for odx, o in enumerate(self.Xo):
                    this_po = p_o[odx]
                    if this_po > 0:
                        this_pio = p_io[self.state_index([i,o])]
                        mi += this_pio*np.log(this_pio / (this_pi * this_po))
        return mi

    def get_all_Cts_R1(self, N, t_max, Nt):
        # right now this is only for comm channels with 2 subsystems
        pvals = np.linspace(0,1,N+1)
        sumto1 = lambda v: sum(v) == 1
        all_pdists_i = np.round(np.array(list(filter(sumto1,list(product(pvals, repeat=self.num_istates))))), 5)
        # print(all_pdists_i)
        num_pdists_i = all_pdists_i.shape[0]

        num_ts = Nt+1
        ts = np.linspace(0,t_max,num_ts)
        tstep = ts[1]

        all_Is_forall_pis = np.zeros((num_pdists_i, num_ts))

        for pind, p_i in tqdm(enumerate(all_pdists_i)):
            all_ps = np.zeros((num_ts, self.num_states))
            all_pis = np.zeros((num_ts, self.num_istates))
            all_pos = np.zeros((num_ts, self.num_ostates))
            all_Is = np.zeros(num_ts)

            # set init p to have output as delta function at self.Xo[0]
            p_init = np.zeros(self.num_states)
            for ival in self.Xi:
                ind = np.intersect1d(np.where(self.X[:,self.i] == ival)[0], np.where(self.X[:,self.o] == self.Xo[0])[0])
                prob = p_i[self.sub_state_index(self.i, ival)]
                np.put(p_init, ind, prob)

            p = p_init
            p_o = self.get_marginal([self.o], self.Xo , self.X, p)
            all_ps[0] = p
            all_pis[0] = p_i
            all_pos[0] = p_o
            all_Is[0] = self.I(p, p_i, p_o)
            for i,t in enumerate(ts[1:]):
                p = expm_multiply((self.Ksys.tocsc()*tstep).tocsr(), p)
                pi = self.get_marginal([self.i], self.Xi , self.X, p)
                po = self.get_marginal([self.o], self.Xo , self.X, p)
                all_ps[i+1] = p
                all_pis[i+1] = pi
                all_pos[i+1] = po
                all_Is[i+1] = self.I(p, pi, po)

            all_Is_forall_pis[pind] = all_Is

        inds = np.argmax(all_Is_forall_pis, axis=0)
        print(np.amax(all_Is_forall_pis, axis=0))
        print(np.argmax(all_Is_forall_pis, axis=0))
        print(all_pdists_i[inds])

    def step (self, p, tstep, K=None):
        p_bef = p
        if K == None:
            p = expm_multiply((self.Ksys.tocsc()*tstep).tocsr(), p)
        else:
            p = expm_multiply((K.tocsc()*tstep).tocsr(), p)
        p_aft = p
        L1 = sum([np.abs(x-y) for x,y in zip(p_bef, p_aft)])
        return L1, p_aft

    def get_Cness_R1(self, N, timestep, L1_thresh, thismodelname):
        # right now this is only for comm channels with 2 subsystems
        pvals = np.linspace(0,1,N+1)
        sumto1 = lambda v: sum(v) == 1
        all_pdists_i = np.round(np.array(list(filter(sumto1,list(product(pvals, repeat=self.num_istates))))), 5)
        num_pdists_i = all_pdists_i.shape[0]

        all_pdists = np.zeros((num_pdists_i, self.num_states))
        all_Is = np.zeros(num_pdists_i)

        all_Is_dict = {}
        for pind, p_i in tqdm(enumerate(all_pdists_i)):

            # set init p to have output as delta function at self.Xo[0]
            p_init = np.zeros(self.num_states)
            for ival in self.Xi:
                ind = np.intersect1d(np.where(self.X[:,self.i] == ival)[0], np.where(self.X[:,self.o] == self.Xo[0])[0])
                prob = p_i[self.sub_state_index(self.i, ival)]
                np.put(p_init, ind, prob)

            p = p_init
            p_o = self.get_marginal([self.o], self.Xo , self.X, p)
            all_pdists[pind] = p

            L1 = 10000
            while L1 > L1_thresh:
                L1, p = step(p, timestep)

            pi = self.get_marginal([self.i], self.Xi , self.X, p)
            po = self.get_marginal([self.o], self.Xo , self.X, p)

            Ival = self.I(p, pi, po)
            all_Is[pind] = Ival
            all_Is_dict[tuple(np.round(N*p_i, 3))] = Ival

        ind = np.argmax(all_Is)
        val = np.amax(all_Is)

        if self.num_istates == 3:
            scale = N
            fig, ax = plt.subplots()
            figure, tax = ternary.figure(ax=ax, scale=scale)
            tax.heatmap(all_Is_dict, vmin = np.amin(all_Is), vmax=val, style = 'hexagonal', cmap = 'viridis')
            tax.boundary(linewidth=2.0)
            tax.set_title(r'$C^{NESS}$')
            # tax.show()
            plt.savefig('{}-Cness.png'.format(thismodelname))
            plt.close()

        self.Cness = val

        return all_pdists[ind], val

    def get_Ki_from_diag(self, ks, L=None, inds=None):
        if L is None:
            L = self.num_istates
            inds = [np.where(self.X[:,0] == l)[0] for l in range(L)]
        kdiag = np.zeros(self.num_states)
        for i in range(L):
            np.put(kdiag, inds[i], ks[i])
        iflips = self.where_spin_i_flips(self.i)
        kdiag_sparse = diags(kdiag)
        K = kdiag_sparse.dot(iflips.dot(kdiag_sparse))
        K = K + diags(np.array(-1*K.sum(axis = 0))[0])
        return K

    def get_Cness_NR(self, N, timestep, L1_thresh):
        L = self.num_istates
        inds = [np.where(self.X[:,0] == l)[0] for l in range(L)]
        all_ks = np.array(list(product(*[np.arange(N+1) for l in range(L)])))
        all_ks = np.delete(all_ks, 0, axis = 0)
        num_ks = all_ks.shape[0]

        all_Is = np.zeros(num_ks)
        all_pinits = np.zeros((num_ks, self.num_states))
        all_pis = np.zeros((num_ks, self.num_istates))
        all_Ks = []

        for n, ks in tqdm(enumerate(all_ks)):

            K = self.get_Ki_from_diag(ks, L, inds)
            init_dist = np.zeros(shape=(self.num_states, ))
            np.put(init_dist, np.where(self.X[:,self.o]==0)[0], 1.0)
            p = init_dist / np.sum(init_dist)

            L1 = 10000
            while L1 > L1_thresh:
                L1, p = self.step(p, timestep, K)

            pi = self.get_marginal([self.i], self.Xi , self.X, p)
            pinit = np.zeros(shape=(self.num_states, ))
            np.put(pinit, np.where(self.X[:,self.o]==0)[0], pi)

            self.set_mechs_RMs(K_input=K)
            self.set_init_dist('test', pinit)

            self.get_to_NESS(timestep = timestep, L1_thresh = L1_thresh)
            all_Is[n] = self.I()
            all_pis[n] = pi
            all_pinits[n] = pinit
            all_Ks.append(K)

        print(all_Is)
        ind = np.argmax(all_Is)
        val = np.amax(all_Is)
        print()
        print(val)
        print(all_pis[ind])
        print(all_Ks[ind].toarray())

        # print(np.round(all_pis, 3))
        print(all_pis[all_Is > 0.9*val])
        print(all_ks[ind])

        self.Cness = val

        return all_Ks[ind], all_pinits[ind]


class Model_R1 (Comm_Channel):

    def set_mechs_RMs (self):

        #input doesn't change state after being set initially here
        mechs = []
        mechs = self.append_output_mechs(mechs)
        self.set_mechs(mechs)
        self.set_all_RMs()
        self.invalid_state_inds = np.array([], dtype=int)


class Model_NR (Comm_Channel):

    def set_mechs_RMs (self, K_input=None):

        if K_input is None:
            mechs = [(self.i, self.sr)]
        else:
            mechs = [(self.i, self.sr*K_input)]
        mechs = self.append_output_mechs(mechs)
        self.set_mechs(mechs)
        self.set_all_RMs()
        self.invalid_state_inds = np.array([], dtype=int)


def generate_all_trajs(Model, modelname, num_runs, num_steps, ness_runtime, param, param_vals, energy_fcns, To_vals, ids, num_msgs):
    '''
        blah

    '''

    def generate_trajs(energy_fcn, To_val, id, num_msg, param=None, param_val=None):

        spc = np.arange(num_msg)
        spin_spaces = {0: spc, 1: spc}
        attrs = {   'H': energy_fcn,
                    'To': [1, To_val]}

        if param is not None:
            this_model_name = '{}_{}-{}_{}-{}_{}-{}_{}-{}_{}-{}'.format(modelname, 'H', energy_fcn, 'To', To_val, param, param_val, 'id', id, 'L', num_msg)
            attrs[param] = param_val
        else:
            this_model_name = '{}_{}-{}_{}-{}_{}-{}_{}-{}_{}-{}'.format(modelname, 'H', energy_fcn, 'To', To_val, 'sr', 0.0, 'id', id, 'L', num_msg)
            attrs['sr'] = 0.0

        print(this_model_name)
        model = Model((1,2), spin_spaces, toroidal = False, input_spin = 0, output_spin = 1, energy_func = energy_fcn, ness_runtime = ness_runtime)
        model.set_attrs(attrs)

        '''
        From tests and knowledge of how energy functions determine
        the input distribution that gives the channel capacity
        '''
        if energy_fcn == 'delta':
            ks = np.ones(model.num_istates)
            init_dist = np.zeros(shape=(model.num_states, ))
            np.put(init_dist, np.where(model.X[:,model.o]==0)[0], 1.0)
            Cness_pinit = init_dist / np.sum(init_dist)

        elif energy_fcn == 'abs' or energy_fcn == 'harmonic':
            ks = np.zeros(model.num_istates)
            ks[0] = 1
            ks[-1] = 1
            whereout0 = np.where(model.X[:,model.o]==0)[0]
            whereinfirst= np.intersect1d(np.where(model.X[:,model.i]==model.Xi[0])[0], whereout0)
            whereinlast= np.intersect1d(np.where(model.X[:,model.i]==model.Xi[-1])[0], whereout0)
            init_dist = np.zeros(shape=(model.num_states, ))
            np.put(init_dist, whereinfirst, 1.0)
            np.put(init_dist, whereinlast, 1.0)
            Cness_pinit = init_dist / np.sum(init_dist)

        Cness_Ki = model.get_Ki_from_diag(ks)
        if modelname == 'model_nr':
            model.set_mechs_RMs(K_input=Cness_Ki)
        elif modelname == 'model_r1':
            model.set_mechs_RMs()

        model.set_init_dist('Cness', Cness_pinit)
        model.get_to_NESS(timestep = 0.01, L1_thresh = 10**(-10))
        model.Cness = model.I()

        '''
        For determining the input distribution that gives the channel capacity
        explicitly, use the following commented out code
        # if modelname == 'model_nr':
        #     Cness_Ki, Cness_pinit = model.get_Cness_NR(N=5, timestep = 1, L1_thresh = 10**(-10))
        #     model.set_mechs_RMs(K_input=attrs['sr']*Cness_Ki)
        #     model.set_init_dist('Cness', Cness_pinit)
        #     model.get_to_NESS(timestep = 0.01, L1_thresh = 10**(-10))
        #     model.Cness = model.I()
        #
        # elif modelname == 'model_r1':
        #     Cness_pinit, Cness = model.get_Cness_R1(N=32, timestep = 1, L1_thresh = 10**(-10), thismodelname=this_model_name)
        #     print(Cness_pinit, Cness)
        #     model.set_init_dist(id, Cness_pinit)
        #     model.get_to_NESS(timestep = 0.01, L1_thresh = 10**(-10))
            # model.get_all_Cts(N=15,t_max=1,Nt=100)
        '''

        # print(np.round(model.Ksys.toarray(),2))
        # print(model.p_init)
        # print(model.p_ness)

        model.write_data(this_model_name, num_runs = num_runs, NESS = True)
        if modelname != 'model_r1':
            model.write_data(this_model_name, num_runs = num_runs, wind = True)

        model.write_data(this_model_name, num_steps = num_steps, pdists = True)

        pf = tb.open_file(this_model_name + '_pdists_relax.h5', mode='r')
        pfr = pf.root
        ts = pfr.pdists_relax.pdists_relax.col('time')
        ps = pfr.pdists_relax.pdists_relax.col('p')[:,0]
        model.write_data(this_model_name, num_runs = num_runs, relax = True, saved_ts = ts, saved_ps = ps)
        pf.close()

        # model.sanity_check(this_model_name)

        print()

    if param is not None:
        for param_val in param_vals:
            for To_val, id, num_msg, energy_fcn in product(To_vals, ids, num_msgs, energy_fcns):
                generate_trajs(energy_fcn, To_val, id, num_msg, param, param_val)

    else:
        for To_val, id, num_msg, energy_fcn in product(To_vals, ids, num_msgs, energy_fcns):
            generate_trajs(energy_fcn, To_val, id, num_msg)


################################################################################
# ALL THE RUNS

energy_fcns = ['delta', 'abs']
To_vals = np.round(np.linspace(0.5, 5.5, 3), 2)
sr_vals = np.round(np.linspace(0.5, 5.5, 3), 2)
initdists = ['Cness']
n_msgs = np.linspace(2, 5, 2).astype(int)

n_runs = 10
n_steps = 200
ness = 10

np.set_printoptions(suppress=True)

# generate_all_trajs(Model_R1, 'model_r1', n_runs, n_steps, ness, None, [], energy_fcns, To_vals, initdists, n_msgs)
# generate_all_trajs(Model_NR, 'model_nr', n_runs, n_steps, ness, 'sr', sr_vals, energy_fcns, To_vals, initdists, n_msgs)

''' red
check paper notes "in the morning" to see how to proceed
remember - #1 priority is to have something down for APS, and right now that's
"what is the behavior of basic comm channels running at their channel capacity?"
+ "show that TUR for info flow holds at NESS"
---> even if it's not very satisfying to you rn, this isn't your endgame,
---> plus it will get a lot of people interested since people are obsessed with how efficient
        channels are at the channel capacity, e.g. as evidenced by SNR vs channel capacity plots
        red in order to PLOT, you have to
            red red THINK
            of what plots you can make to mimic what you're seeing in literature
            for example, in addition to plots you already outlined, you can do
                (To and/or sr) vs C as dots, where color of dots gives value of a given quantity like EP
                    accrued during relaxation or for a fixed time in the NESS
            red red red READ
            Thomas and Cover on channel capacity and general intuitions and reasons to look into channel capacity
            (remember you have to prepare for possible questions!)
            for example, think of channel capacity and power efficiency limit

---> also it will be important for you to have started reading and digesting
        papers for TD of distributed computational systems before APS
        so that you can have meaningful discussions with people there

and only after that whole talk is 100% done and prepared,
and whether or not you'll end up putting it into APS talk,
you can look into other interesting things like, what happens to thermodynamic behavior of the channel
when it's run below its channel capacity
'''

################################################################################
# ONE RUN FOR TESTING

# spin_spaces = {0: [0,1,2], 1: [0,1,2]}
# attrs = {}
# attrs['xi'] = 1
# attrs['mu'] = 1
# attrs['eta'] = 1
# attrs['eps'] = 1
# attrs['sr'] = 1
# attrs['spo'] = 1
# attrs['To'] = [1, 1]
# attrs['Ti'] = [1]
#
# print('test')
# test = Model_NR((1,2), spin_spaces, toroidal = False, input_spin = 0, output_spin = 1, ness_runtime = 10)
# test.set_attrs(attrs)
# test.set_mechs_RMs()
# print(np.round(test.Ksys.toarray(),2))
# test.set_init_dist('uni')
# test.get_to_NESS(timestep = 0.01, L1_thresh = 10**(-10))
# print(test.p_init)
# print(test.p_ness)
#
# test.write_data('test', num_runs = 10, NESS = True, wind = True)
# test.write_data('test', num_steps = 1000, pdists = True)
# pf = tb.open_file('test' + '_pdists_relax.h5', mode='r')
# pfr = pf.root
# ts = pfr.pdists_relax.pdists_relax.col('time')
# ps = pfr.pdists_relax.pdists_relax.col('p')[:,0]
# test.write_data('test', num_runs = 10, relax = True, saved_ts = ts, saved_ps = ps)
# pf.close()

# test.sanity_check('test')
