import os
import tables as tb
import numpy as np
from numpy import savez_compressed as save_array
from scipy.sparse import save_npz

from copy import deepcopy

import networkx as nx
from scipy.stats import norm, binom, randint, uniform
from scipy.stats import bernoulli, betabinom, binom, boltzmann, dlaplace, geom, hypergeom, logser, nbinom, planck, poisson, randint, skellam, yulesimon, zipf, zipfian
from scipy.sparse.linalg import expm, expm_multiply, eigs
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix, csc_matrix
from scipy.sparse import identity, diags, kron
from itertools import combinations, product, permutations, pairwise
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

'''
red
Need to fix the ragged array creation in the rate matrix construction code, just to get rid of the error.
'''
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


'''
**********************************************************************************************************************************
*** Note that a mechanism is identified through the set of coordinates in which it causes simultaneous state transitions -
    This doesn't have to map a priori to a reservoir, even though it could. This is relevant, e.g., for coarse graining.

Goals:
We want to be able to simulate stochastic processes as composite processes so as to be able to track / make use of
    Quantities
        Information flows between puppets of different mechanisms
            Looking at the overlap of puppets will also be interesting
            E.g., ATP may in a biochemical network be a hub node of the hypergraph
        EP due to each mechanism
        EF due to each mechanism
            Question: when can we use the EF instead of the EP to bound quantities (e.g., to be experimentally accessible)
        KL Divergences D(p_Lv || pi_Lv) where pi_Lv is the steady state distribution of the leaders of v
            need to build into this code the ability to obtain the steady state distributions for any subset of coordinates
            to do that, recognize that you can find the steady state distribution for the smallest unit in which L(v) is contained
            and then take the marginal
        Conservative and non-conservative forces
    Descriptions of structure
        Dependency graph incidence matrices (if relevant)
            both the |N| x |N| formulation and the |N| x |E| formulation
            their Laplacians
            eigenvalues of their Laplacians,
            and associated quantities like, spectral gaps, etc.
        Hypergraph-mechanism (rule) incidence matrices
            their Laplacians
            eigenvalues of their Laplacians,
            and associated quantities like, spectral gaps, etc.
            red
            note that for hypergraphs you can define the notion of N-back leaders, which is the concept of leaders but extended further
                to look at leaders-of-leaders, or leaders-of-leaders-of-leaders, and so on... This might have some significance, not sure.
        Units and !!! can define incidence matrix over unit structures!
            their Laplacians
            eigenvalues of their Laplacians,
            and associated quantities like, spectral gaps, etc.
        Rate matrix graph incidence matrices
            their Laplacians
            and in particular ratio of imaginary vs real parts of the eigenvalues
            (of course have to set a particular method to set the RMs,
            and check that it's not due to that method but due to the structure of the DAG / hypergraph)
    Results
        Units (within each unit and/or involving multiple units)
            mix-and-match TURs [of information flows]
            speed limit theorems
            decomposition of information flows
            decompositions of EPs
            relations to KL divergences
that apply to them.

Plans for how to generate data and interface with it:
I need to check how much data Pytables should be requiring for the trajectories and make sure it's not being inflated.
**********************************************************************************************************************************
'''

'''
**********************************************************************************************************************************
Setting up some non-class functions
**********************************************************************************************************************************
'''

def translate(a, d):
    '''Convert all items in array from a --> b using a dictionary map d'''
    return np.vectorize(d.__getitem__)(a)

def list_of_lists_translate(l, d):
    '''Convert all items in list of lists from a --> b using a dictionary map d'''
    return [[d[a] for a in sublist] for sublist in l]

def dict_translate(d_to_translate, d):
    '''Convert all values of dict from a --> b using a dictionary map d'''
    return {key: translate(value, d) for key, value in d_to_translate.items()}

def union(arr):
    '''Return multi-set union of a set of sets'''
    out = arr[0]
    for l in arr[1:]:
        out = out + l
    return np.array(out)

def overlap(l1, l2):
    '''Return intersection of two sets'''
    return list(set(l1) & set(l2))

def intersection (arr):
    '''Return intersection of a set of sets'''
    l = len(arr)
    out = arr[0]
    if l > 1:
        out = overlap(arr[0], arr[1])
        for i in range(l):
            if i > 1:
                out = overlap(out, arr[i])
    return out

def union_minus_intersection(arr):
    '''Return union minus intersection of a set of sets'''
    u = union(arr)
    i = intersection(arr)
    return [x for x in u if x not in i]

rng = np.random.default_rng()

'''
**********************************************************************************************************************************
Setting up named tuples and classes that store sets of variables that will be referenced repeatedly. Instances of named tuples are
immutable, but instances of classes are, of course, mutable.
**********************************************************************************************************************************
'''
ValidNexts = namedtuple('ValidNexts', ('indices, states, propensities'))
Indices = namedtuple('Indices', ('ind', 'inds_i', 'inds_pups', 'inds_leads', 'inds_negleads', 'inds_units'))
Quant_Dots = namedtuple('Quant_Dots', ('q_dot', 'if_io_dot', 'if_oi_dot', 'mi_dot', 'mi_inst', 'zdot'))

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


'''
**********************************************************************************************************************************
**********************************************************************************************************************************
Main constructor class for composite processes
**********************************************************************************************************************************
**********************************************************************************************************************************
'''

class Composite_Process():

    '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Given that we now understand that CPs are at the end of the day defined by rules, we need to find a way to translate
    the hypergraph representing the graph structure of the rules governing the allowed system transitions to network types
    like the small world, etc....

    We need to understand

    (1) How an MPP's dependency graph maps onto its directed multihypergraph of rules
        (there will be an ensemble of DMHGs consistent with any given dependency graph).
        --> It might be the case that you can convert the DMHG over nodes to be such that you have
            bidirectional links between any two of the puppets of a reservoir, and
            unidirectional links wherever there are catalysts --
            can there be unidirectional links in other cases?

    (2) If the CP is not an MPP, then what would it mean for it to be a
        small world / stochasticblock / modular / hierarchical
        - it's going to relate to some property of the incidence matrix of the DMHG
        (and not the adjacency matrix over complexes because that erases information about the underlying nodes)
        - e.g., if a generic CP is modular, then it will be hyperedges that give the same property:
            higher in-cluster connectivity than between-cluster connectivity

    (3) Do you have to write your own code to generate CPs of a given type (small world, moduler, hierarchical, etc.?
        Or can you somehow modify the networkx stuff to translate to hypergraph???? That probably depends how well your
        random sampling reflects the underlying ?degree distribution?)
        Maybe use networkx to generate the DG of style (1) and then randomly convert that into a DMHG given
        some other parameters, like n_i to sample from, etc.

        --> This is going to be a fairly subtle and complex task so come back to this later.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    '''

    def set_coordinates (self, coords):
        '''
        inputs:
            coords (arr): of coordinates, can be strings or just an array that is the list of its indices
        '''
        self.coordinates = coords
        self.indices_from_coords = {coord: index for index, coord in enumerate(coords)}
        self.coords_from_indices = {index: coord for index, coord in enumerate(coords)}
        self.num_coords = len(coords)
        self.coordinates_indices = np.arange(self.num_coords)

    def set_spaces (self, X_i):
        '''
        inputs:
            X_i (dict): {i: state space of coordinate i}
        creates:
            X_u (dict): {u: state space of unit u}
            X (arr of arr): of all possible joint states of system
        '''
        self.X_i = X_i
        self.X_i_dims = {i: len(xi) for i, xi in X_i.items()}
        self.X_p = {v: self.get_joint_state_space({i: self.X_i[i] for i in pups_v}) for v, pups_v in self.puppets.items()}
        self.X_l = {v: self.get_joint_state_space({i: self.X_i[i] for i in leaders_v}) for v, leaders_v in self.leaders.items()}
        self.X_nl = {v: self.get_joint_state_space({i: self.X_i[i] for i in negleaders_v}) for v, negleaders_v in self.negleaders.items()}
        self.X_u = {uind: self.get_joint_state_space({i: self.X_i[i] for i in u}) for uind, u in self.units.items()}
        self.X = self.get_joint_state_space(self.X_i)
        self.num_states = len(self.X)

    def get_joint_state_space (self, coord_spaces):
        '''
        inputs:
            coord_spaces: (dict) {i: (arr) of coordinate i's state space for all desired coordinates}
        outputs:
            (arr) of all states in joint state space
                Note that the order in which product() iterates through states is what the later efficient rate matrix
                construction algorithm relies on in order to generate the correct rate matrices
        '''
        return np.array([np.array(s) for s in product(*coord_spaces.values())])

    def set_mechanisms (self, mechs):
        '''
        inputs:
            mechs (arr): of mechanisms, can be strings or just an array that is the list of its indices, i.e., np.arange(|V|)
        '''
        self.mechanisms = mechs
        self.num_mechs = len(mechs)
        # self.indices_from_mechs = {mech: index for index, mech in enumerate(mechs)}
        # self.mechs_from_indices = {index: mech for index, mech in enumerate(mechs)}

    def set_puppets (self, puppets):
        '''
        inputs:
            puppets (dict): { v : np.array([i for i in P(v)]) }
        '''
        self.puppets = puppets
        # self.puppets_indices = dict_translate(self.puppets, self.indices_from_coords)

    def set_leaders (self, leaders):
        '''
        inputs:
            leaders (dict): { v : np.array([i for i in L(v)]) }
        '''
        self.leaders = leaders
        # self.leaders_indices = dict_translate(self.leaders, self.indices_from_coords)
        self.negleaders = {v: np.array([i for i in self.coordinates if i not in vleads]) for v, vleads in self.leaders.items()}
        # self.negleaders_indices = dict_translate(self.negleaders, self.indices_from_coords)

    def set_rules (self, rules):
        '''
            rules (dict): { v : np.array([[n_i for i in r] for r in R(v)]) }
        '''
        self.rules = rules

    '''
    For now, just set the units yourself.
    red red red
    Later, create an algorithm that will obtain units assuming that there is a
    tight unit structure, as David calls it in his paper.
    '''

    def set_units (self, units):
        '''
            units (dict): { u : np.array([i for i in U(u)]) }
        '''
        self.units = units
        self.num_units = len(self.units)

    def get_dmhg_incidence(self):
        '''
        What would this look like? |N| x |R|

            A_ij = n_i(j) if i in P(j), for node i and rule j
            Note that since every forward transition rule will have a corresponding backward transition rule with
                n_i* = - n_i
            this A matrix will have linearly dependent columns, so when extracting the eigenvalues/vectors, just
            choose only the forward transitions?

        The adjacency matrix over complexes would also make sense to include because then you could actually put the
        rates in it... Again, come back to this. Not needed at this precise moment (6/6/2022)
        '''
        pass

    def construct_all_Mvs (self):
        '''
        creates:
            Ms_v (list of scipy sparse matrices): each for a given mechanism v,
                with ones where the puppets of v change state, and zeros otherwise
        '''

        def reduced_matrix_constructor(space_i, n_i):
            dim = len(space_i)
            out = np.zeros((dim, dim))
            for i1, a in enumerate(space_i):
                for i2, b in enumerate(space_i):
                    if a - b == n_i:
                        out[i1][i2] = 1
            return out

        def iterative_hypergraph_constructor(coord_spaces, puppets, rules):
            m = []
            for rule in rules:
                _m_ = []
                for i, X_i in reversed(list(coord_spaces.items())):
                    dim_i = len(X_i)
                    if i in puppets:
                        which_puppet = np.where(puppets == i)[0][0]
                        _m_.append(coo_matrix(reduced_matrix_constructor(X_i, rule[which_puppet])))
                    else:
                        _m_.append(identity(dim_i))
                m_ = kron(_m_[1], _m_[0])
                for m1 in _m_[2:]:
                    m_ = kron(m1, m_)
                m.append(m_)

            '''If m can be written as a sparse matrix then the following can be a single line sum over axis 0'''
            if len(m) == 1:
                return m[0]
            elif len(m) == 2:
                return m[0] + m[1]
            else:
                out = m[0] + m[1]
                for mat in m[2:]:
                    out = out + mat
                return out

        self.Ms_v = [iterative_hypergraph_constructor(self.X_i, self.puppets[mech], self.rules[mech]) for mech in self.mechanisms]

        def check_symmetric(a, rtol=1e-05, atol=1e-05):
            return np.allclose(a, a.T, rtol=rtol, atol=atol)

        print ('This should be an array of True only', [check_symmetric(M.toarray()) for M in self.Ms_v])

    def construct_rate_diags_fcn_of_state (self, fcn, v):
        '''
        inputs:
            fcn (fcn): that takes as input the state, the mechanism, and a dictionary of the leaders of all mechanisms,
                and returns a float that is placed along a diagonal matrix, the entry corresponding to the given state
        creates:
            (scipy sparse diags matrix): with fcn(state) for every state in the joint system state space
        '''
        return diags(np.array([fcn(x, v, self.leaders) for x in self.X]))

    def construct_all_RMs (self, bef_state_fcns, aft_state_fcns):
        '''
        inputs:
            bef_state_fcn (fcn): of type fcn in self.construct_rate_diags_fcn_of_state, for the before state of the transition
            aft_state_fcn (fcn): of type fcn in self.construct_rate_diags_fcn_of_state, for the after state of the transition
        creates:
            Ks_v: (list of scipy sparse arrs) containing the K(v) for all mechanisms v of the system
            Ksys: (scipy sparse arr) representing the system rate matrix (= sum_v K(v))
            all of these have size |self.X| x |self.X|
        '''
        self.construct_all_Mvs()
        Ks_v = []
        c = 0
        for v_ind, v in enumerate(self.mechanisms):
            transitions_v = self.Ms_v[v_ind]

            if isinstance(bef_state_fcns, list):
                diags_v_bef_stubs = [self.construct_rate_diags_fcn_of_state(bef_fcn, v) for bef_fcn in bef_state_fcns]
                diags_v_aft_stubs = [self.construct_rate_diags_fcn_of_state(aft_fcn, v) for aft_fcn in aft_state_fcns]
                K_v_stubs = [diags_v_bef.dot(transitions_v.dot(diags_v_aft)) for diags_v_bef, diags_v_aft in zip(diags_v_bef_stubs, diags_v_aft_stubs)]
                # print(self.Ms_v)
                K_v = diags(np.zeros(self.num_states))
                for K_ in K_v_stubs:
                    K_v += K_
                K_v = K_v + diags(np.array(-1*K_v.sum(axis = 0))[0])
                Ks_v.append(K_v)
            else:
                diags_v_bef = self.construct_rate_diags_fcn_of_state(bef_state_fcns)
                diags_v_aft = self.construct_rate_diags_fcn_of_state(aft_state_fcns)

                K_v = diags_v_bef.dot(transitions_v.dot(diags_v_aft))
                K_v = K_v + diags(np.array(-1*K_v.sum(axis = 0))[0])
                Ks_v.append(K_v)

            if c == 0:
                K = K_v
            else:
                K = K + K_v

            c += 1

        self.Ks_v = Ks_v
        self.Ksys = K

        Ks_v_dense = [K_v.toarray() for K_v in Ks_v]

    def coord_state_index (self, coord, coord_state):
        '''
        inputs:
            sub: (int) indicating coordinate in self.sub_state_spaces
            coord_state: state of coordinate in question
        outputs:
            (int) indicating index of coord_state in sub's state space
        '''
        return np.where(self.X_i[coord]==coord_state)[0][0]

    def state_index (self, state):
        '''
        inputs:
            state: (arr) representing state of joint system
        outputs:
            (int) indicating index of state in joint state space
        '''
        return reduce(np.intersect1d, [np.where(self.X[:,i]==s)[0] for i,s in enumerate(state)])[0]

    def state_indices(self, mstate, msubset, jstate_spc):
        '''
        inputs:
            mstate: (arr) representing state of a subset of all the coordinates involved in jstate_spc
            msubset: (arr) representing which coordinates, out of the ones involved in jstate_spc, are in the subset
            jstate_spc: (arr of arr) containing all joint states from which we want to pluck. At a minimum it contains the
                coordinates in msubset.
        outputs:
            (arr) of indices in jstate_spc corresponding to where the msubset takes on the state mstate
        '''
        if len(msubset) > 1:
            return reduce(np.intersect1d, [np.where(jstate_spc[:,sub]==mstate[i])[0] for i,sub in enumerate(msubset)])
        else:
            return reduce(np.intersect1d, [np.where(jstate_spc[:,sub]==mstate)[0] for i,sub in enumerate(msubset)])

    def subset_state_index(self, mstate, mstate_spc):
        return self.state_indices(mstate, np.arange(len(mstate)), mstate_spc)[0]

    def get_marginal(self, subset, marg_spc, joint_spc, joint_pdist):
        '''
        inputs:
            subset: (arr) representing which coordinates, out of the ones involved in joint_spc, are in the subset
            marg_spc: (arr of arr) of all joint states of subset
            joint_spc: (arr of arr) containing all joint states from which we want to pluck. At a minimum it contains the
                coordinates in subset
            joint_pdist: (arr) joint probability distribution of joint_spc
        outputs:
            (arr) marginal probability distribution of the subset of coordinates
        '''
        return np.array([np.sum(np.take(joint_pdist, self.state_indices(ms, subset, joint_spc))) for ms in marg_spc])

    def set_init_dist(self, id, p):
        '''
        sets system properties:
            self.id (str) indicating the type of the initial distribution
            self.p_init (arr) initial probability distribution to use when generating trajectories
        '''
        self.id = id
        self.p_init = p

        '''
        red red red
        probably would be good to create some methods of setting the initial distribution according to some general types of
        distributions, e.g., uniform, gaussian, etc., as you did in the cell ligand-receptor-memory model code
        '''

    def sample_init_dist(self, init_dist = None):
        '''
        inputs:
            init_dist: (arr) provide this iff you want to use something that is not the system's initial distribution, for example
                as may be desired to explore the simplex for the distribution that gives the channel capacity at the NESS
        outputs:
            index: (int) indicated index of randomly sampled state
            state: (arr) indicating the value of the sampled state
        '''
        if init_dist is None:
            index = np.random.choice(self.num_states, p=self.p_init)
        else:
            index = np.random.choice(self.num_states, p=init_dist)
        state = self.X[index]
        return index, state

    def get_to_NESS (self, timestep, L1_thresh, ness_runtime):
        '''
        Iterates pnew = p*e^(self.Ksys*timestep) repeatedly until L1(pnew, p) < L1_thresh

        inputs:
            timestep: (float)
            L1_thresh: (float)

        sets system properties:
            self.time_to_NESS: (float) time required to reach the NESS
            self.p_ness: (arr) probability distribution of the joint system at the NESS
            self.ps_i_ness: (dict of sub: arr) marginal probability distribution of each sub at the NESS
            self.p_io_ness: (arr) marginal probability distribution of the joint state of the input and output at the NESS
        '''
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
        self.ps_i_ness = {i: self.get_marginal([i_ind], self.X_i[i], self.X, self.p_ness) for i_ind, i in enumerate(self.coordinates)}
        self.ness_runtime = ness_runtime

        return t, self.p_ness, self.ps_i_ness


    def get_state_tau_ago (self, traj, curr_t, tau):
        target = curr_t - tau
        return traj['states'][np.where(traj['times'] < target)[0][-1]]


    def get_state_inds (self, state):
        '''
        Returns indices of
            the joint state in the joint state space
        '''
        joint_si = self.state_index(state)
        coord_si = [self.coord_state_index(self.coordinates[coord], state[coord]) for coord in self.coordinates_indices]
        pupp_si = [self.subset_state_index(state[pups], self.X_p[v]) for v, pups in self.puppets.items()]
        leader_si = [self.subset_state_index(state[leads], self.X_l[v]) for v, leads in self.leaders.items()]
        negleader_si = [self.subset_state_index(state[negleads], self.X_nl[v]) for v, negleads in self.negleaders.items()]
        unit_si = [self.subset_state_index(state[unit], self.X_u[uind]) for uind, unit in self.units.items()]
        inds = Indices(joint_si, coord_si, pupp_si, leader_si, negleader_si, unit_si)
        return inds

    def get_options (self, ind):
        '''
        Obtains options for next state when generating trajectories via the Gillespie algorithm

        inputs:
            ind: (int) index of current joint state

        outputs:
            all_propensities: (arr) of the propensities [rate matrix elements] of switching to every possible other state,
                due to each possible (i,v) pair
            valid_nexts: (ValidNexts) of (arr, arr, arr) of indices of valid next states, valid next states, and the propensities
                of transitioning to those next states (valid states are any states for which the propensity to switch to it from
                the current state is greater than 0)
            sum_propensities: (float) sum of the propensities of all valid next states
        '''
        all_propensities = np.concatenate([K_v.getcol(ind).toarray().T[0] for K_v in self.Ks_v], axis=0)
        valid_next_indices = np.where(all_propensities > 0)
        valid_next_states = self.X[np.mod(valid_next_indices[0], self.num_states)]
        valid_propensities = all_propensities[all_propensities > 0]
        sum_propensities = np.sum(valid_propensities)
        valid_nexts = ValidNexts(valid_next_indices, valid_next_states, valid_propensities)
        # print(valid_nexts)
        return all_propensities, valid_nexts, sum_propensities

    def get_sojourn(self, sum_propensities):
        '''
        Returns a random sojourn time based on the total propensity to leave the current state
        '''
        r = rng.random()
        return -1 * np.log(r) / sum_propensities

    def get_transition(self, valid_nexts, sum_propensities):
        '''
        Based on the current state, the valid next states, and the total propensity to leave the current state, implements the
        Gillespie algorithm for randomly choosing the next state
        '''
        valid_next_indices, valid_next_states, valid_propensities = valid_nexts
        propensity_scale = np.cumsum(valid_propensities)/sum_propensities
        r = rng.random()
        choice = np.where(propensity_scale > r)[0].min()
        next_state = valid_next_states[choice]
        next_state_ind = self.state_index(next_state)
        which_v = valid_next_indices[0][choice] // self.num_states
        return next_state_ind, next_state, which_v


    def set_traj_inits(self, p_init=None):
        '''
        For cleaner code, separates off all initializations required for generating a trajectory
        '''
        if p_init is None:
            p = self.p_init
        else:
            p = p_init
        init_index, init_state = self.sample_init_dist(p)
        init_indices = self.get_state_inds(init_state)
        t = 0
        curr_state_ind = init_index
        next_state_ind = None

        return t, init_indices, curr_state_ind, next_state_ind, p


    def row_append_pdist_relax(self, row, time, pdist):
        '''
        Appends a row to a PyTables table tracking the evolution of the joint and marginal probability distributions at evenly
        spaced intervals throughout the relaxation process
        '''
        row['time'] = time
        row['pdist'] = pdist
        row.append()

    def row_append_ness(self, row, *args):
        '''
        Appends a row to a PyTables table for every state transition occurring in a trajectory, for an NESS
        '''
        if len(args) == 8:
            time, ind_v, state_ind, state_i_inds, state_pup_inds, state_lead_inds, state_neglead_inds, state_u_inds = args
        elif len(args) == 3:
            time, ind_v, state_inds = args
            state_ind, state_i_inds, state_pup_inds, state_lead_inds, state_neglead_inds, state_u_inds = state_inds

        row['time'] = time
        row['which_v'] = ind_v
        row['state_ind'] = state_ind
        row['state_i_inds'] = state_i_inds
        row['state_pup_inds'] = state_pup_inds
        row['state_lead_inds'] = state_lead_inds
        row['state_neglead_inds'] = state_neglead_inds
        row['state_u_inds'] = state_u_inds

        row.append()

    def row_append_rlx(self, row, *args):
        '''
        Appends a row to a PyTables table for every state transition occurring in a trajectory, for relaxation
        '''
        if len(args) == 9:
            time, ind_v, state_ind, state_i_inds, state_pup_inds, state_lead_inds, state_neglead_inds, state_u_inds, pdist = args
        elif len(args) == 4:
            time, ind_v, state_inds, pdist = args
            state_ind, state_i_inds, state_pup_inds, state_lead_inds, state_neglead_inds, state_u_inds = state_inds

        row['time'] = time
        row['which_v'] = ind_v
        row['state_ind'] = state_ind
        row['state_i_inds'] = state_i_inds
        row['state_lead_inds'] = state_lead_inds
        row['state_neglead_inds'] = state_neglead_inds
        row['state_u_inds'] = state_u_inds
        row['pdist'] = pdist

        row.append()

    '''
    ------------------------------------------------------------------------------------------------------------------------------
    Below functions are for generating trajectories or probability distributions during NESS or relaxation processes
    ------------------------------------------------------------------------------------------------------------------------------
    '''

    def gentraj_NESS (self, row):
        '''
        Given the row iterator for the PyTables file in which the trajectory will be stored,
        generate one NESS trajectory and store it
        '''
        t, init_indices, curr_state_ind, next_state_ind, p = self.set_traj_inits(self.p_ness)
        self.row_append_ness(row, t, -2, init_indices)

        while t < self.ness_runtime:
            if t > 0:
                curr_state_ind = next_state_ind

            all_propensities, valid_nexts, sum_propensities = self.get_options(curr_state_ind)
            if sum_propensities == 0:
                print()
                print('wtf! and the time is ', t)
                print(self.X[curr_state_ind])
                print(all_propensities)
                print(self.get_sojourn(sum_propensities))
            sojourn = self.get_sojourn(sum_propensities)
            t += sojourn

            if t < self.ness_runtime:
                next_state_ind, next_state, which_v = self.get_transition(valid_nexts, sum_propensities)
                next_inds = self.get_state_inds(next_state)
                self.row_append_ness(row, t, which_v, next_inds)

            else:
                inds = self.get_state_inds(self.X[curr_state_ind])
                self.row_append_ness(row, self.ness_runtime, -2, inds)


    def generate_pdists_relaxation(self, row, num_steps):
        '''
        Given the row iterator for the PyTables file in which the evolution of probability distribution during relaxation
        will be stored, generate those probability distributions and store them
        '''
        p = self.p_init
        t = 0
        self.row_append_pdist_relax(row, t, p)
        timestep = self.time_to_NESS / num_steps

        def step(p, t):
            p_bef = p
            p = expm_multiply((self.Ksys.tocsc()*timestep).tocsr(), p)
            p_aft = p
            L1 = sum([np.abs(x-y) for x,y in zip(p_bef, p_aft)])
            t += timestep
            return L1, p_aft, t

        while t <= self.time_to_NESS:
            L1, p, t = step(p, t)
            self.row_append_pdist_relax(row, t, p)

    def gentraj_relax (self, row, saved_ts, saved_ps):
        '''
        Given the row iterator for the PyTables file in which the trajectory will be stored,
        generate one relaxation process trajectory and store it
        '''
        t, init_indices, curr_state_ind, next_state_ind, p = self.set_traj_inits()
        self.row_append_rlx(row, t, -2, init_indices, p)

        while t < self.time_to_NESS:
            if t > 0:
                curr_state_ind = next_state_ind

            all_propensities, valid_nexts, sum_propensities = self.get_options(curr_state_ind)
            sojourn = self.get_sojourn(sum_propensities)
            t += sojourn

            if t < self.time_to_NESS:
                nearest_t_ind = np.where(saved_ts <= t)[0][-1]
                extra_bitof_t = t - saved_ts[nearest_t_ind]
                p = saved_ps[nearest_t_ind]
                # print(saved_ps)
                # print(p)
                # print()
                # print(self.X)
                p = expm_multiply((self.Ksys.tocsc()*extra_bitof_t).tocsr(), p)
                # print(p)
                # print()
                # print(self.X_l)
                next_state_ind, next_state, which_v = self.get_transition(valid_nexts, sum_propensities)
                next_inds = self.get_state_inds(next_state)
                self.row_append_rlx(row, t, which_v, next_inds, p)

            else:
                nearest_t_ind = np.where(saved_ts <= self.time_to_NESS)[0][-1]
                extra_bitof_t = self.time_to_NESS - saved_ts[nearest_t_ind]
                p = saved_ps[nearest_t_ind]
                p = expm_multiply((self.Ksys.tocsc()*extra_bitof_t).tocsr(), p)
                inds = self.get_state_inds(self.X[curr_state_ind])
                self.row_append_rlx(row, self.time_to_NESS, -2, inds, p)

    '''
    ------------------------------------------------------------------------------------------------------------------------------
    Below functions are for more PyTables logistics and making it easy to batch-generate data of trajectories, distributions, etc.
    ------------------------------------------------------------------------------------------------------------------------------
    '''
    def save_attributes (self):
        '''
        Save important attributes that will be needed to go from raw trajectories to thermodynamic calculations
        '''
        np.save('coordinates.npy', self.coordinates)
        np.save('mechanisms.npy', self.mechanisms)

        np.save('puppets.npy', self.puppets)
        np.save('leaders.npy', self.leaders)
        np.save('negleaders.npy', self.negleaders)
        np.save('units.npy', self.units)

        np.save('X.npy', self.X)
        np.save('X_i.npy', self.X_i)
        np.save('X_u.npy', self.X_u)
        np.save('X_p.npy', self.X_p)
        np.save('X_l.npy', self.X_l)
        np.save('X_nl.npy', self.X_nl)

        np.save('p_ness.npy', self.p_ness)
        np.save('p_init.npy', self.p_init)

        np.save('Ksys.npy', self.Ksys)
        np.save('Ks_v.npy', self.Ks_v)

        try:
            np.save('dg_incidence.npy', self.dg_incidence)
        except:
            pass
        try:
            np.save('dg_adjacency.npy', self.dg_adjacency)
        except:
            pass
        try:
            np.save('hg_incidence.npy', self.dg_incidence)
        except:
            pass
        try:
            np.save('coord_leaders.npy', self.coord_leaders)
        except:
            pass

    def set_group (self, f, root, typ):
        '''
        Set the attributes to store in the group of a given PyTables (specified by the file f and its root)
            with typ (str) indicating your name of whatever is being stored in the group.
            These group attributes indicate the parameters that were used to generate the
                [trajectories, distributions, etc.] that are being stored in the group
        '''
        group = f.create_group(root, typ)
        return group

    def set_cols(self):
        '''
        Set the different columns for the different types of PyTables files that could be generated
        '''

        class cols_ness(tb.IsDescription):
            time = tb.Float64Col(pos=0)
            which_v = tb.Float64Col(pos=1)
            state_ind = tb.Int16Col(pos=2)
            state_i_inds = tb.Int16Col(shape=(self.num_coords, ), pos=3)
            state_pup_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=4)
            state_lead_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=5)
            state_neglead_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=6)
            state_u_inds = tb.Int16Col(shape=(self.num_units, ), pos=7)

        class cols_rlx(tb.IsDescription):
            time = tb.Float64Col(pos=0)
            which_v = tb.Float64Col(pos=1)
            state_ind = tb.Int16Col(pos=2)
            state_i_inds = tb.Int16Col(shape=(self.num_coords, ), pos=3)
            state_pup_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=4)
            state_lead_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=5)
            state_neglead_inds = tb.Int16Col(shape=(self.num_mechs, ), pos=6)
            state_u_inds = tb.Int16Col(shape=(self.num_units, ), pos=7)
            pdist = tb.Float64Col(shape = (self.num_states, ), pos=8)

        class cols_pdist(tb.IsDescription):
            time = tb.Float64Col(pos=0)
            pdist = tb.Float64Col(shape = (self.num_states, ), pos=1)


        self.cols_ness = cols_ness
        self.cols_rlx = cols_rlx
        self.cols_pdist = cols_pdist


    def write_pdists(self, fname, num_steps):
        '''
        inputs:
            fname: (str) name of PyTables file being used to store probability distributions for the current set of conditions
                the rows of this file store the probability distributions (joint and several marginal ones) for each timestep
            cols: (obj) class of type tb.IsDescription that defines the columns for the table data
            num_steps: (int) number of timesteps for which the probability distributions will be calculated and stored

        output is that by the end of the function running, the file with fname will contain all the desired data
        '''
        tb.file._open_files.close_all()
        f = tb.open_file(fname, mode='w')
        root = f.root
        grp = self.set_group(f, root, 'pdists_rlx')
        table = f.create_table(grp, 'pdists_rlx', self.cols_pdist, 'pdists_rlx')
        r = table.row
        self.generate_pdists_relaxation(r, num_steps)
        table.flush()
        f.close()

    def write_traj(self, fname, NorR, cols, num_runs, saved_ts=None, saved_ps=None):
        '''
        inputs:
            fname: (str) name of PyTables file being used to store trajectories for the current set of conditions
                the rows of this file store the probability distributions (joint and several marginal ones) for each timestep
            NorR: (str) eithr 'n' or 'r' to indicate either NESS or relaxation process to be run
            cols: (obj) class of type tb.IsDescription that defines the columns for the trajectory data
            num_runs: (int) number of trajectories to generate
            if relaxation trajectory:
            saved_ts: (arr) of timesteps that were stored in the corresponding pdists_relax file for this set of conditions
            saved_ps: (arr of arr) of probability distributions that were stored in the corresponding pdists_relax file for this
                set of conditions

        output is that by the end of the function running, the file with fname will contain all the desired simulation data
        '''
        def run(runfile, root, group, run_fcn, *args):
            for i in trange(num_runs):
                this_run = 'run_' + str(i).zfill(7)
                table = runfile.create_table(group, this_run, cols, this_run)
                r = table.row
                run_fcn(r, *args)
                table.flush()

        f = tb.open_file(fname, mode='w')
        root = f.root
        grp = self.set_group(f, root, 'runs')
        if NorR == 'n':
            run(f, root, grp, self.gentraj_NESS)
        elif NorR == 'r':
            run(f, root, grp, self.gentraj_relax, saved_ts, saved_ps)
        f.close()


    def write_raw_data(self, c, dirname, num_runs = 0, num_steps = 0, NESS = False, relax = False, pdists = False, saved_ts = None, saved_ps = None):
        '''
        Wrapper function used to generate and write any of the possible types of trajectories, probability distributions, etc.
        inputs already described in the functions that require them
        '''
        self.set_cols()

        if os.path.isdir(dirname):
            os.chdir(dirname)
        else:
            os.mkdir(dirname)
            os.chdir(dirname)

        self.save_attributes()

        if NESS:
            print('Writing NESS trajectories...')
            self.write_traj(c + '_traj_NESS.h5', 'n', self.cols_ness, num_runs)
        if relax:
            print('Writing relaxation trajectories...')
            self.write_traj(c + '_traj_rlx.h5', 'r', self.cols_rlx, num_runs, saved_ts, saved_ps)
        if pdists:
            print('Writing probability distributions during relaxation path...')
            self.write_pdists(c + '_pdists_rlx.h5', num_steps)

        os.chdir('..')


'''
Example instantiation of a composite process
'''
cp = Composite_Process()
cp.set_coordinates(np.arange(4))
cp.set_mechanisms(np.arange(3))
test_pups = {0: np.array([0,2]), 1: np.array([1,2]), 2: np.array([1,3])}
cp.set_puppets(test_pups)
cp.set_leaders(test_pups)
# test_rules = {0: np.array([[1, -1, 2], [-1, 1, -2]]),
#               1: np.array([[2, -1], [-2, 1]]),
#               2: np.array([[2, -2], [-2, 2]])}
test_rules = {i: np.array([[1, -1], [-1, 1]]) for i in cp.mechanisms}
cp.set_rules(test_rules)
test_units = {0: np.array([0,1,2,3]), 1: np.array([0])}
cp.set_units(test_units)
cp.set_spaces({0: np.array([0,1,2]), 1: np.array([0,1,2,3,4]), 2: np.array([0,1,2,3]), 3: np.array([0,1,2])})
cp.construct_all_Mvs()


def match_filter(x, v, leaders):
    '''
    If min of an arr is equal to max of the arr then the arr contains only one unique element

    inputs:
        x (arr): state
        v (str or int): mechanism
        leaders (dict): {v (str or int): L(v) (arr)}
    '''
    return int( np.min(x[leaders[v]]) == np.max(x[leaders[v]]) )

def mismatch_filter(x, v, leaders):
    '''
    inputs:
        x (arr): state
        v (str or int): mechanism
        leaders (dict): {v (str or int): L(v) (arr)}
    '''
    return int( np.min(x[leaders[v]]) != np.max(x[leaders[v]]) )

epsilon = 0.3

'''
Currently this depends on the type of transition (matching or mismatching) and
the only dependence on the leaders are to get the marginal state to check that they match or mismatch after
'''

'''
Here the g_{x_{L(v)}} is simply np.exp( +- 0.5*len(np.unique(x[leaders[v]])))
Only need one fcn here because we want anyway for the rate into/out of a mismatched state out of/into a matched one
to be symmetric and dependent on the number of mismatched coordinates in the mismatched state
'''
def mismatch_rates(x, v, leaders):
    return (mismatch_filter(x, v, leaders))*np.exp(-0.5*len(np.unique(x[leaders[v]])))

'''
And for mismatched state --> different mismatched state the rate is just 1.
Which is the geometric mean of the match --> mismatch and mismatch --> match transitions
'''

'''
red red red
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TO DO:

(1) Create the calculate_cp_quantities functions, essentially converting get_dots and write_quants to this new paradigm
    --> make sure to build in quants that can be used to test the speed limits and TURs for information flows

(2) Test with catalysts (what does it actually look like when implementing?)
    --> for example when the catalyst is a regulated conductor, then the voltage of that conductor will change the "amount"
    of charge being "conserved", see paper notes in LA CENTRAL notebook
(3) Calculate unit-based quantities; which means also defining units (we will assume tight unit structures)
(3) Set up an example of N = 3 pbits (in detail!) and test
(4) Set up an example of a simple CRN with ATP and test
(5) Modify plotting functions so you can plot things like info flows and speed limits and TURs for info flows etc.

(6) Fill in the functions for creating the appropriate incidence / adjacency matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
red red red
'''

cp.construct_all_RMs([match_filter, mismatch_rates, mismatch_filter], [mismatch_rates, match_filter, mismatch_filter])

# print(cp.Ksys)
# print()

this_model_name = 'test_ckt'
this_dirname = this_model_name

N = 4
states_with_charge_N = np.array([1 if sum(x) == N else 0 for x in cp.X])
cp.set_init_dist('uniform_chargeN', states_with_charge_N/sum(states_with_charge_N))
cp.get_to_NESS(timestep = 1, L1_thresh = 10**-10, ness_runtime = 50)
cp.write_raw_data(this_model_name, this_dirname, num_runs = 50, NESS = True)
cp.write_raw_data(this_model_name, this_dirname, num_steps = 100, pdists = True)

pf = tb.open_file(this_dirname + '\\' + this_model_name + '_pdists_rlx.h5', mode='r')
pfr = pf.root
ts = pfr.pdists_rlx.pdists_rlx.col('time')
ps = pfr.pdists_rlx.pdists_rlx.col('pdist')

cp.write_raw_data(this_model_name, this_dirname, num_runs = 50, relax = True, saved_ts = ts, saved_ps = ps)
pf.close()

'''
red red red
Later, write wrapper functions that generate trajectories for all kinds of different systems
'''



'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Think about this issue:

Note that chemical reaction networks are going to be a weird case. They are represented as directed MULTIhypergraphs.
OR as directed graphs over the complexes.

So, the underlying DAG is... ? nonexistent? Or perhaps it is the pairwise DAG that can be uniquely associated with the directed
hypergraph? So how would the concept of the leaders of a mechanism (chemical reaction) be? I think:
    vertices are species
    puppets of the chemical reaction v are the species involved. sources are the reactants, targets are the products
        net/effective chem rxns modelled so that you don't have species that are both source and target
    leaders of the chemical reaction v are either:
        (1) all other species involved in chem rxns with puppets of v, OR
        (2) all other species which are sources (reactants) in chem rxns where a puppet of v is a target (product)
        - probably first one because the concentrations of all species involved in a chem reaction affect the free energy of rxn

So for the chem rxn you would actually start with the hypergraph and then produce the dependency graph
    (where all dependencies are bidirectional between i and j if i and j reside in a hyperedge together?)

Okay, now thinking of a circuit, as in Freitas & Esposito, you don't really have an underlying DAG - so what is the concept of
leaders in this case I think hyperneighbors of hyperneighbors... (checking David's paper) Doesn't actually specify it but this is what would
make sense for this context.

(3) Make functions that get the eigenvalues of the M(v) matrices, and see how it affects final Ksys matrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
