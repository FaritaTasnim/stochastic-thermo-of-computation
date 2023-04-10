import os
import numpy as np
import tables as tb
from tqdm import trange, tqdm

from collections import namedtuple

from scipy.stats import norm, binom, randint, uniform
from scipy.stats import bernoulli, betabinom, binom, boltzmann, dlaplace, geom, hypergeom, logser, nbinom, planck, poisson, randint, skellam, yulesimon, zipf, zipfian
from scipy.sparse.linalg import expm, expm_multiply, eigs
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix, csc_matrix
from scipy.sparse import identity, diags, kron
from itertools import combinations, product, permutations, pairwise
from functools import reduce


def choose_two (n):
    return int(0.5*n*(n-1))

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Have to grab lots of stuff, e.g., num_coords, num_mechs, and all the npy files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

System_Properties = namedtuple('System_Properties', (   'coordinates',
                                                        'num_coords',
                                                        'mechanisms',
                                                        'num_mechs',
                                                        'puppets',
                                                        'leaders',
                                                        'negleaders',
                                                        'units',
                                                        'num_units',
                                                        'X',
                                                        'X_i',
                                                        'X_u',
                                                        'X_p',
                                                        'X_l',
                                                        'X_nl',
                                                        'p_ness',
                                                        'p_init',
                                                        'Ksys',
                                                        'Ks_v'))

Iter_Params = namedtuple('Iter_Params', (   'coords_expanded',
                                            'X_i_expanded',
                                            'coordsP2',
                                            'coordsC2',
                                            'all_X_cP2',
                                            'all_X_cC2',
                                            'mechsP2',
                                            'mechsC2',
                                            'puppetsP2',
                                            'leadersP2',
                                            'puppetsC2',
                                            'leadersC2',
                                            'puppetsP2_U',
                                            'leadersP2_U',
                                            'puppetsC2_U',
                                            'leadersC2_U',
                                            'all_X_pP2',
                                            'all_X_lP2',
                                            'all_X_pC2',
                                            'all_X_lC2',
                                            'unit_indsP2',
                                            'unit_indsC2',
                                            'unitsP2',
                                            'unitsC2',
                                            'unitsP2_U',
                                            'unitsC2_U',
                                            'all_X_uP2',
                                            'all_X_uC2'))

Iter_Dists = namedtuple('Iter_Dists', ( 'coord_dists',
                                        'coordsP2_dists',
                                        'coordsC2_dists',
                                        'puppet_dists',
                                        'leader_dists',
                                        'negleader_dists',
                                        'puppetsP2_dists',
                                        'leadersP2_dists',
                                        'puppetsC2_dists',
                                        'leadersC2_dists',
                                        'unit_dists',
                                        'unitsP2_dists',
                                        'unitsC2_dists'))

State_Indices = namedtuple('State_Indices', (   'bef_state_ind',
                                                'bef_state_i_inds',
                                                'bef_state_pup_inds',
                                                'bef_state_lead_inds',
                                                'bef_state_neglead_inds',
                                                'bef_state_u_inds',
                                                'aft_state_ind',
                                                'aft_state_i_inds',
                                                'aft_state_pup_inds',
                                                'aft_state_lead_inds',
                                                'aft_state_neglead_inds',
                                                'aft_state_u_inds',
                                                'bef_state_cP2_inds',
                                                'bef_state_cC2_inds',
                                                'bef_state_pP2_inds',
                                                'bef_state_lP2_inds',
                                                'bef_state_pC2_inds',
                                                'bef_state_lC2_inds',
                                                'bef_state_uP2_inds',
                                                'bef_state_uC2_inds',
                                                'aft_state_cP2_inds',
                                                'aft_state_cC2_inds',
                                                'aft_state_pP2_inds',
                                                'aft_state_lP2_inds',
                                                'aft_state_pC2_inds',
                                                'aft_state_lC2_inds',
                                                'aft_state_uP2_inds',
                                                'aft_state_uC2_inds'))


'''
red red red
Also make a variable class for quantities ?
'''

all_quant_names = [ 'time',
                    #
                    'if_mechs',
                    'if_coords',
                    'if_puppets',
                    'if_leaders',
                    'if_units',
                    #
                    'mi_inst_mechs',
                    'mi_inst_coords',
                    'mi_inst_puppets',
                    'mi_inst_leaders',
                    'mi_inst_units',
                    #
                    'mi_tot_mechs',
                    'mi_tot_coords',
                    'mi_tot_puppets',
                    'mi_tot_leaders',
                    'mi_tot_units',
                    #
                    'ep',
                    'zeta_mechs',
                    'ep_units',
                    #
                    'ef',
                    'ef_mechs',
                    'ef_puppets',
                    'ef_leaders',
                    'ef_units',
                    #
                    's',
                    's_coords',
                    's_puppets',
                    's_leaders',
                    's_units']

def quants_row_append (r, quants):
    '''
    Let quants be an array of values for all of the quantities in cols_qs for a given timestep/entry
    '''
    for qname, q in zip(all_quant_names, quants):
        r[qname] = q
    r.append()

def get_state_indices(mstate, mcoords, joint_spc):
    '''
    inputs:
        mstate: (arr) representing state of a subset of all the coordinates involved in jstate_spc
        mcoords: (arr) representing which coordinates, out of the ones involved in jstate_spc, are in the subset
        jstate_spc: (arr of arr) containing all joint states from which we want to pluck. At a minimum it contains the
            coordinates in mcoords.
    outputs:
        (arr) of indices in jstate_spc corresponding to where the msubset takes on the state mstate
    '''
    try:
        return reduce(np.intersect1d, [np.where(joint_spc[:,sub]==mstate[i])[0] for i,sub in enumerate(mcoords)])
    except:
        return reduce(np.intersect1d, np.where(joint_spc[:,mcoords]==mstate)[0])

def get_joint_state_space (coord_spaces):
    '''
    inputs:
        coord_spaces: (dict) {i: (arr) of coordinate i's state space for all desired coordinates}
    outputs:
        (arr) of all states in joint state space
            Note that the order in which product() iterates through states is what the later efficient rate matrix
            construction algorithm relies on in order to generate the correct rate matrices
    '''
    return np.array([np.array(s) for s in product(*coord_spaces.values())])

def get_marginal (marg_coords, marg_spc, joint_spc, joint_dist):
    '''
    inputs:
        marg_coords: (arr) representing which coordinates, out of the ones involved in joint_spc, are in the subset
        marg_spc: (arr of arr) of all joint states of marg_coords
        joint_spc: (arr of arr) containing all joint states from which we want to pluck. At a minimum it contains the
            coordinates in marg_coords
        joint_pdist: (arr) joint probability distribution of joint_spc
    outputs:
        (arr) marginal probability distribution of the subset of subsystems
    '''
    return np.array([np.sum(np.take(joint_dist, get_state_indices(marg_state, marg_coords, joint_spc))) for marg_state in marg_spc])

'''The below quantity calculation functions are for NESS only -- need separate ones for relaxation'''

def get_setA_minus_setB (A, B):
    return np.array([i for i in A if i not in B])

def get_if_dot (pAB, pB, aftAB, befAB, aftB, befB):
    '''
    Calculates the trajectory-level change in I^{coords_A --> coords_B} due to state transition (bef --> aft)
    '''
    if_dot = np.log((pAB[aftAB]*pB[befB])/(pAB[befAB]*pB[aftB]))

    return if_dot

def get_mi_inst (pAB, pA, pB, stateAB, stateA, stateB):
    '''
    Calculates the trajectory-level mutual information I(coords_A; coords_B) for a given state and joint distribution
    '''

    mi_inst = np.log(pAB[stateAB]/(pA[stateA]*pB[stateB]))

    return mi_inst

def get_ef_dot (bef, aft, which_v, Ks_v):
    '''
    Calculates the entropy flow from the relevant reservoir/mechanism into the system due to state transition (bef -->(v) aft)
    '''
    if which_v == -2:
        ef_dot = 0
    else:
        ef_dot = np.log(Ks_v[which_v][aft,bef]/Ks_v[which_v][bef,aft])
    return ef_dot

def get_s (pcoords, statecoords):
    '''
    Calculates the trajectory-level entropy of a given state from a given joint distribution (or marginal if coords != system)
    '''
    return -1*np.log(pcoords[statecoords])

def get_zeta_dot (bef, aft, bef_state_Lv, aft_state_Lv, joint_dist, which_v, Ks_v, pLv):
    '''
    Calculates the EP due to a particular transition (this is NESS so you don't have to account for the [nonexistent] change of the dist)
    '''
    ef_dot = get_ef_dot(bef, aft, which_v, Ks_v)
    delta_s = get_s(pLv, aft_state_Lv) - get_s(pLv, bef_state_Lv)
    if_dot = get_if_dot(joint_dist, pLv, aft, bef, aft_state_Lv, bef_state_Lv)
    ep_dot = ef_dot + delta_s + if_dot
    return ef_dot, delta_s, if_dot, ep_dot

def get_P2 (pool):
    return [[a, b] for a, b in permutations(pool, 2)]

def get_C2 (pool):
    return [[a, b] for a, b in combinations(pool, 2)]

def get_unions_of_pairs (arr_of_pairs):
    return [np.unique(np.concatenate(pair)) for pair in arr_of_pairs]

def get_X_subsets (X_i, subsets):
    return [get_joint_state_space({i: X_i[i] for i in subset}) for subset in subsets]

def get_dists_subsets (subsets, X_subsets, X, joint_dist):
    return [get_marginal(subset, X_subset, X, joint_dist) for subset, X_subset in zip(subsets, X_subsets)]

def get_state_inds (X, jstate_ind, subsets, X_subsets):
    return [get_state_indices(X[jstate_ind][subset], np.arange(len(subset)), X_subset)[0] for subset, X_subset in zip(subsets, X_subsets)]

def read_run_row (row):
    time = row['time']
    state_ind = row['state_ind']
    state_i_inds = row['state_i_inds']
    state_pup_inds = row['state_pup_inds']
    state_lead_inds = row['state_lead_inds']
    state_neglead_inds = row['state_neglead_inds']
    state_u_inds = row['state_u_inds']
    return time, state_ind, state_i_inds, state_pup_inds, state_lead_inds, state_neglead_inds, state_u_inds

def get_all_state_inds (X, jstate_ind, ips):
    state_cP2_inds = get_state_inds(X, jstate_ind, ips.coordsP2, ips.all_X_cP2)
    state_cC2_inds = get_state_inds(X, jstate_ind, ips.coordsC2, ips.all_X_cC2)
    state_pP2_inds = get_state_inds(X, jstate_ind, ips.puppetsP2_U, ips.all_X_pP2)
    state_lP2_inds = get_state_inds(X, jstate_ind, ips.leadersP2_U, ips.all_X_lP2)
    state_pC2_inds = get_state_inds(X, jstate_ind, ips.puppetsC2_U, ips.all_X_pC2)
    state_lC2_inds = get_state_inds(X, jstate_ind, ips.leadersC2_U, ips.all_X_lC2)
    state_uP2_inds = get_state_inds(X, jstate_ind, ips.unitsP2_U, ips.all_X_uP2)
    state_uC2_inds = get_state_inds(X, jstate_ind, ips.unitsC2_U, ips.all_X_uC2)
    return state_cP2_inds, state_cC2_inds, state_pP2_inds, state_lP2_inds, state_pC2_inds, state_lC2_inds, state_uP2_inds, state_uC2_inds

def calc_NESS_qs (write_table, run, Ks, ips, ids, sp):
    '''
    Probably best to loop through the run once and at every transition, step the tabulation of each of the desired quantities
    in the same way that I had been doing get_dots in the previous methodology
    '''

    def get_if_dots (si):
        if_coords_dot = np.array([get_if_dot(dist_cP2, ids.coord_dists[cP2[1]], aft_state_cP2, bef_state_cP2, si.aft_state_i_inds[cP2[1]], si.bef_state_i_inds[cP2[1]]) for cP2, dist_cP2, bef_state_cP2, aft_state_cP2 in zip(ips.coordsP2, ids.coordsP2_dists, si.bef_state_cP2_inds, si.aft_state_cP2_inds)])
        if_puppets_dot = np.array([get_if_dot(dist_pP2, ids.puppet_dists[vP2[1]], aft_state_pP2, bef_state_pP2, si.aft_state_pup_inds[vP2[1]], si.bef_state_pup_inds[vP2[1]]) for vP2, dist_pP2, bef_state_pP2, aft_state_pP2 in zip(ips.mechsP2, ids.puppetsP2_dists, si.bef_state_pP2_inds, si.aft_state_pP2_inds)])
        if_leaders_dot = np.array([get_if_dot(dist_lP2, ids.leader_dists[vP2[1]], aft_state_lP2, bef_state_lP2, si.aft_state_lead_inds[vP2[1]], si.bef_state_lead_inds[vP2[1]]) for vP2, dist_lP2, bef_state_lP2, aft_state_lP2 in zip(ips.mechsP2, ids.leadersP2_dists, si.bef_state_lP2_inds, si.aft_state_lP2_inds)])
        if_units_dot = np.array([get_if_dot(dist_uP2, ids.unit_dists[uP2[1]], aft_state_uP2, bef_state_uP2, si.aft_state_u_inds[uP2[1]], si.bef_state_u_inds[uP2[1]]) for uP2, dist_uP2, bef_state_uP2, aft_state_uP2 in zip(ips.unit_indsP2, ids.unitsP2_dists, si.bef_state_uP2_inds, si.aft_state_uP2_inds)])
        return if_coords_dot, if_puppets_dot, if_leaders_dot, if_units_dot

    def get_mi_insts (si):
        mi_inst_mechs = np.array([get_mi_inst(sp.p_ness, dist_Lv, dist_Lv_, si.bef_state_ind, si.bef_state_lead_inds[mech], si.bef_state_neglead_inds[mech]) for mech, (dist_Lv, dist_Lv_) in enumerate(zip(ids.leader_dists, ids.negleader_dists))])
        mi_inst_coords = np.array([get_mi_inst(dist_cC2, ids.coord_dists[cC2[0]], ids.coord_dists[cC2[1]], bef_state_cC2, si.bef_state_i_inds[cC2[0]], si.bef_state_i_inds[cC2[1]]) for cC2, dist_cC2, bef_state_cC2 in zip(ips.coordsC2, ids.coordsC2_dists, si.bef_state_cC2_inds)])
        mi_inst_puppets = np.array([get_mi_inst(dist_pC2, ids.puppet_dists[vC2[0]], ids.puppet_dists[vC2[1]], bef_state_pC2, si.bef_state_pup_inds[vC2[0]], si.bef_state_pup_inds[vC2[1]]) for vC2, dist_pC2, bef_state_pC2 in zip(ips.mechsC2, ids.puppetsC2_dists, si.bef_state_pC2_inds)])
        mi_inst_leaders = np.array([get_mi_inst(dist_lC2, ids.leader_dists[vC2[0]], ids.leader_dists[vC2[1]], bef_state_lC2, si.bef_state_lead_inds[vC2[0]], si.bef_state_lead_inds[vC2[1]]) for vC2, dist_lC2, bef_state_lC2 in zip(ips.mechsC2, ids.leadersC2_dists, si.bef_state_lC2_inds)])
        mi_inst_units = np.array([get_mi_inst(dist_uC2, ids.unit_dists[uC2[0]], ids.unit_dists[uC2[1]], bef_state_uC2, si.bef_state_u_inds[uC2[0]], si.bef_state_u_inds[uC2[1]]) for uC2, dist_uC2, bef_state_uC2 in zip(ips.unit_indsC2, ids.unitsC2_dists, si.bef_state_uC2_inds)])
        return mi_inst_mechs, mi_inst_coords, mi_inst_puppets, mi_inst_leaders, mi_inst_units

    def get_delta_ss (si):
        s_coords_bef = np.array([get_s(cdist, cstate) for cdist, cstate in zip(ids.coord_dists, si.bef_state_i_inds)])
        s_puppets_bef = np.array([get_s(pupdist, pupstate) for pupdist, pupstate in zip(ids.puppet_dists, si.bef_state_pup_inds)])
        s_leaders_bef = np.array([get_s(leaddist, leadstate) for leaddist, leadstate in zip(ids.leader_dists, si.bef_state_lead_inds)])
        s_units_bef = np.array([get_s(udist, ustate) for udist, ustate in zip(ids.unit_dists, si.bef_state_u_inds)])

        s_coords_aft = np.array([get_s(cdist, cstate) for cdist, cstate in zip(ids.coord_dists, si.aft_state_i_inds)])
        s_puppets_aft = np.array([get_s(pupdist, pupstate) for pupdist, pupstate in zip(ids.puppet_dists, si.aft_state_pup_inds)])
        s_leaders_aft = np.array([get_s(leaddist, leadstate) for leaddist, leadstate in zip(ids.leader_dists, si.aft_state_lead_inds)])
        s_units_aft = np.array([get_s(udist, ustate) for udist, ustate in zip(ids.unit_dists, si.aft_state_u_inds)])

        delta_s_coords = s_coords_bef - s_coords_aft
        delta_s_puppets = s_puppets_bef - s_puppets_aft
        delta_s_leaders = s_leaders_bef - s_leaders_aft
        delta_s_units = s_units_bef - s_units_aft

        return delta_s_coords, delta_s_puppets, delta_s_leaders, delta_s_units

    def get_ef_dots (ef_dot, v):
        ef_mechs_dot = np.array([ef_dot if v != -2 and v_ == v else 0 for v_ in sp.mechanisms])
        ef_puppets_dot = np.array([ef_dot if v != -2 and pv in sp.puppets[v] else 0 for pv in sp.puppets.values()])
        ef_leaders_dot = np.array([ef_dot if v != -2 and lv in sp.puppets[v] else 0 for lv in sp.leaders.values()])
        ef_units_dot = np.array([ef_dot if v != -2 and u in sp.puppets[v] else 0 for u in sp.units.values()])
        return ef_mechs_dot, ef_puppets_dot, ef_leaders_dot, ef_units_dot

    '''Initializing all quantities to zero'''
    num_coordsC2 = choose_two(sp.num_coords)
    num_mechsC2 = choose_two(sp.num_mechs)
    num_unitsC2 = choose_two(sp.num_units)
    quants = np.array([ 0,
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (2*num_coordsC2, )),
                        np.zeros(shape = (2*num_mechsC2, )),
                        np.zeros(shape = (2*num_mechsC2, )),
                        np.zeros(shape = (2*num_unitsC2, )),
                        #
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (num_coordsC2, )),
                        np.zeros(shape = (num_mechsC2, )),
                        np.zeros(shape = (num_mechsC2, )),
                        np.zeros(shape = (num_unitsC2, )),
                        #
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (num_coordsC2, )),
                        np.zeros(shape = (num_mechsC2, )),
                        np.zeros(shape = (num_mechsC2, )),
                        np.zeros(shape = (num_unitsC2, )),
                        #
                        0,
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_units, )),
                        #
                        0,
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_units, )),
                        #
                        0,
                        np.zeros(shape = (sp.num_coords, )),
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_mechs, )),
                        np.zeros(shape = (sp.num_units, ))
                        ], dtype=object)

    run_len = len(run)
    for b_ind, a_ind in pairwise(np.arange(run_len)):
        bef_row = run[b_ind]
        aft_row = run[a_ind]

        bef_time, bef_state_ind, bef_state_i_inds, bef_state_pup_inds, bef_state_lead_inds, bef_state_neglead_inds, bef_state_u_inds = read_run_row (bef_row)
        aft_time, aft_state_ind, aft_state_i_inds, aft_state_pup_inds, aft_state_lead_inds, aft_state_neglead_inds, aft_state_u_inds = read_run_row (aft_row)

        v = int(aft_row['which_v'])
        sojourn = aft_time - bef_time

        bef_state_cP2_inds, bef_state_cC2_inds, bef_state_pP2_inds, bef_state_lP2_inds, bef_state_pC2_inds, bef_state_lC2_inds, bef_state_uP2_inds, bef_state_uC2_inds = get_all_state_inds (sp.X, bef_state_ind, ips)
        aft_state_cP2_inds, aft_state_cC2_inds, aft_state_pP2_inds, aft_state_lP2_inds, aft_state_pC2_inds, aft_state_lC2_inds, aft_state_uP2_inds, aft_state_uC2_inds = get_all_state_inds (sp.X, aft_state_ind, ips)

        this_si = State_Indices(bef_state_ind,
                                bef_state_i_inds,
                                bef_state_pup_inds,
                                bef_state_lead_inds,
                                bef_state_neglead_inds,
                                bef_state_u_inds,
                                aft_state_ind,
                                aft_state_i_inds,
                                aft_state_pup_inds,
                                aft_state_lead_inds,
                                aft_state_neglead_inds,
                                aft_state_u_inds,
                                bef_state_cP2_inds,
                                bef_state_cC2_inds,
                                bef_state_pP2_inds,
                                bef_state_lP2_inds,
                                bef_state_pC2_inds,
                                bef_state_lC2_inds,
                                bef_state_uP2_inds,
                                bef_state_uC2_inds,
                                aft_state_cP2_inds,
                                aft_state_cC2_inds,
                                aft_state_pP2_inds,
                                aft_state_lP2_inds,
                                aft_state_pC2_inds,
                                aft_state_lC2_inds,
                                aft_state_uP2_inds,
                                aft_state_uC2_inds)

        if_coords_dot, if_puppets_dot, if_leaders_dot, if_units_dot = get_if_dots (this_si)
        mi_inst_mechs, mi_inst_coords, mi_inst_puppets, mi_inst_leaders, mi_inst_units = get_mi_insts (this_si)
        mi_tot_mechs = sojourn * mi_inst_mechs
        mi_tot_coords = sojourn * mi_inst_coords
        mi_tot_puppets = sojourn * mi_inst_puppets
        mi_tot_leaders = sojourn * mi_inst_leaders
        mi_tot_units = sojourn * mi_inst_units
        ef_dot, _, if_dot, ep_dot = get_zeta_dot(bef_state_ind, aft_state_ind, bef_state_lead_inds[v], aft_state_lead_inds[v], sp.p_ness, v, Ks, ids.leader_dists[v])
        if_mechs_dot = np.array([if_dot if v_ == v else 0 for v_ in sp.mechanisms])
        ef_mechs_dot, ef_puppets_dot, ef_leaders_dot, ef_units_dot = get_ef_dots(ef_dot, v)
        s = get_s(sp.p_ness, bef_state_ind)
        delta_s_coords, delta_s_puppets, delta_s_leaders, delta_s_units = get_delta_ss(this_si)
        zeta_mechs_dot = ef_mechs_dot + delta_s_leaders + if_mechs_dot
        ep_units_dot = ef_units_dot + delta_s_units

        quants_dots = np.array([sojourn,
                                if_mechs_dot,
                                if_coords_dot,
                                if_puppets_dot,
                                if_leaders_dot,
                                if_units_dot,
                                mi_inst_mechs,
                                mi_inst_coords,
                                mi_inst_puppets,
                                mi_inst_leaders,
                                mi_inst_units,
                                mi_tot_mechs,
                                mi_tot_coords,
                                mi_tot_puppets,
                                mi_tot_leaders,
                                mi_tot_units,
                                ep_dot,
                                zeta_mechs_dot,
                                ep_units_dot,
                                ef_dot,
                                ef_mechs_dot,
                                ef_puppets_dot,
                                ef_leaders_dot,
                                ef_units_dot,
                                s,
                                delta_s_coords,
                                delta_s_puppets,
                                delta_s_leaders,
                                delta_s_units
                                ], dtype=object)


        quants = quants + quants_dots

        wt_row = write_table.row
        quants_row_append(wt_row, quants)

    print(quants[16])
    print()
    return np.array(quants, dtype=object) # size: num_quants

def to_full(a):
    output = np.full((len(a), max(map(len, a))), np.nan)
    for i, row in enumerate(a):
        output[i, :len(row)] = row
    return output

def calc_rlx_ifs ():
    pass

def calc_rlx_efs ():
    pass

def calc_rlx_s ():
    pass

def calc_rlx_mis ():
    '''
    Do this one first because it's the hardest
    '''
    pass # both for mi_inst and mi_tot

'''
EP calculations will just, I think be a sum of the calculations of EF and S and IF
'''

def calc_rlx_qs (write_table, run, Ks, ips, sp, pr):
    '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    START HERE WHEN BACK ---

    Is there any way to get these distributions in a vectorized fashion,
    so that the operations are faster? Maybe if the operations were done
    in the transpose fashion to the current paradigm...
    !!!
    This is important because we need all these distributions at transition times too,
    for each trajectory...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    '''
    timesteps = pr.pdists_rlx.pdists_rlx.col('time')
    pdists_at_times = pr.pdists_rlx.pdists_rlx.col('pdist')
    print(pdists_at_times.shape)
    pdists_at_times_coords = np.array([get_dists_subsets(ips.coords_expanded, ips.X_i_expanded.values(), sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_coordsP2 = np.array([get_dists_subsets(ips.coordsP2, ips.all_X_cP2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_coordsC2 = np.array([get_dists_subsets(ips.coordsC2, ips.all_X_cC2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_puppets = np.array([get_dists_subsets(sp.puppets.values(), sp.X_p.values(), sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_puppetsP2 = np.array([get_dists_subsets(ips.puppetsP2_U, ips.all_X_pP2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_puppetsC2 = np.array([get_dists_subsets(ips.puppetsC2_U, ips.all_X_pC2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_leaders = np.array([get_dists_subsets(sp.leaders.values(), sp.X_l.values(), sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_leadersP2 = np.array([get_dists_subsets(ips.leadersP2_U, ips.all_X_lP2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_leadersC2 = np.array([get_dists_subsets(ips.leadersC2_U, ips.all_X_lC2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_negleaders = np.array([get_dists_subsets(sp.negleaders.values(), sp.X_nl.values(), sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_units = np.array([get_dists_subsets(sp.units.values(), sp.X_u.values(), sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_unitsP2 = np.array([get_dists_subsets(ips.unitsP2_U, ips.all_X_uP2, sp.X, p) for p in pdists_at_times], dtype=object)
    pdists_at_times_unitsC2 = np.array([get_dists_subsets(ips.unitsC2_U, ips.all_X_uC2, sp.X, p) for p in pdists_at_times], dtype=object)


def calc_qs (c, read_footer, write_footer, sp, cols, rlx = False):
    '''
    This function keeps track of trajectory-level quantities for saved NESS trajectories
    '''

    rf = tb.open_file(c + read_footer, mode = 'r')
    wf = tb.open_file(c + write_footer, mode = 'w')
    rfr = rf.root
    wfr = wf.root
    wfg = wf.create_group(wfr, 'quants')

    coords_expanded = np.expand_dims(sp.coordinates, 1)
    X_i_expanded = {i: np.expand_dims(X_i, 1) for i, X_i in sp.X_i.items()}

    coordsP2 = get_P2(sp.coordinates)
    coordsC2 = get_C2(sp.coordinates)

    all_X_cP2 = get_X_subsets(sp.X_i, coordsP2)
    all_X_cC2 = get_X_subsets(sp.X_i, coordsC2)

    mechsP2 = get_P2(sp.mechanisms)
    mechsC2 = get_C2(sp.mechanisms)

    puppetsP2 = get_P2(sp.puppets.values())
    leadersP2 = get_P2(sp.leaders.values())
    puppetsC2 = get_C2(sp.puppets.values())
    leadersC2 = get_C2(sp.leaders.values())

    puppetsP2_U = get_unions_of_pairs(puppetsP2)
    leadersP2_U = get_unions_of_pairs(leadersP2)
    puppetsC2_U = get_unions_of_pairs(puppetsC2)
    leadersC2_U = get_unions_of_pairs(leadersC2)

    all_X_pP2 = get_X_subsets(sp.X_i, puppetsP2_U)
    all_X_lP2 = get_X_subsets(sp.X_i, leadersP2_U)
    all_X_pC2 = get_X_subsets(sp.X_i, puppetsC2_U)
    all_X_lC2 = get_X_subsets(sp.X_i, leadersC2_U)

    unit_indsP2 = get_P2(sp.units.keys())
    unit_indsC2 = get_C2(sp.units.keys())

    unitsP2 = get_P2(sp.units.values())
    unitsC2 = get_C2(sp.units.values())

    unitsP2_U = get_unions_of_pairs(unitsP2)
    unitsC2_U = get_unions_of_pairs(unitsC2)

    all_X_uP2 = get_X_subsets(sp.X_i, unitsP2_U)
    all_X_uC2 = get_X_subsets(sp.X_i, unitsC2_U)

    iter_params = Iter_Params(  coords_expanded,
                                X_i_expanded,
                                coordsP2,
                                coordsC2,
                                all_X_cP2,
                                all_X_cC2,
                                mechsP2,
                                mechsC2,
                                puppetsP2,
                                leadersP2,
                                puppetsC2,
                                leadersC2,
                                puppetsP2_U,
                                leadersP2_U,
                                puppetsC2_U,
                                leadersC2_U,
                                all_X_pP2,
                                all_X_lP2,
                                all_X_pC2,
                                all_X_lC2,
                                unit_indsP2,
                                unit_indsC2,
                                unitsP2,
                                unitsC2,
                                unitsP2_U,
                                unitsC2_U,
                                all_X_uP2,
                                all_X_uC2)

    if not rlx: #if NESS trajectories
        '''
        red red red
        need to pull sp.p_ness, mechanism, puppets, leaders, units,
        '''

        coord_dists = get_dists_subsets(coords_expanded, X_i_expanded.values(), sp.X, sp.p_ness)

        coordsP2_dists = get_dists_subsets(coordsP2, all_X_cP2, sp.X, sp.p_ness)
        coordsC2_dists = get_dists_subsets(coordsC2, all_X_cC2, sp.X, sp.p_ness)

        puppet_dists = get_dists_subsets(sp.puppets.values(), sp.X_p.values(), sp.X, sp.p_ness)
        leader_dists = get_dists_subsets(sp.leaders.values(), sp.X_l.values(), sp.X, sp.p_ness)
        negleader_dists = get_dists_subsets(sp.negleaders.values(), sp.X_nl.values(), sp.X, sp.p_ness)

        puppetsP2_dists = get_dists_subsets(puppetsP2_U, all_X_pP2, sp.X, sp.p_ness)
        leadersP2_dists = get_dists_subsets(leadersP2_U, all_X_lP2, sp.X, sp.p_ness)
        puppetsC2_dists = get_dists_subsets(puppetsC2_U, all_X_pC2, sp.X, sp.p_ness)
        leadersC2_dists = get_dists_subsets(leadersC2_U, all_X_lC2, sp.X, sp.p_ness)

        unit_dists = get_dists_subsets(sp.units.values(), sp.X_u.values(), sp.X, sp.p_ness)

        unitsP2_dists = get_dists_subsets(unitsP2_U, all_X_uP2, sp.X, sp.p_ness)
        unitsC2_dists = get_dists_subsets(unitsC2_U, all_X_uC2, sp.X, sp.p_ness)

        iter_dists = Iter_Dists(coord_dists,
                                coordsP2_dists,
                                coordsC2_dists,
                                puppet_dists,
                                leader_dists,
                                negleader_dists,
                                puppetsP2_dists,
                                leadersP2_dists,
                                puppetsC2_dists,
                                leadersC2_dists,
                                unit_dists,
                                unitsP2_dists,
                                unitsC2_dists)

        i = 0
        for run in tqdm(rfr.runs):

            this_run = 'run_' + str(i).zfill(7)
            wt = wf.create_table(wfg, this_run, cols, this_run)
            calc_NESS_qs(wt, run, sp.Ks_v, iter_params, iter_dists, sp)
            wt.flush()
            i += 1

    else: # if relaxation trajectories

        pf = tb.open_file(c + '_pdists_rlx.h5', mode='r')
        pr = pf.root

        i = 0
        for run in tqdm(rfr.runs):


            '''
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            red red red
            Fill in here
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            '''
            this_run = 'run_' + str(i).zfill(7)
            wt = wf.create_table(wfg, this_run, cols, this_run)

            iter_dists = None

            calc_rlx_qs(wt, run, sp.Ks_v, iter_params, sp, pr)
            wt.flush()
            i += 1

    rf.close()
    wf.close()


# def write_rlx_qs (c):
#     '''
#
#     '''
#     rf = tb.open_file(c + '_traj_rlx.h5', mode = 'r')
#     wf = tb.open_file(c + '_traj_quants_rlx.h5', mode = 'w')
#     rfr = rf.root
#     wfr = wf.root


def write_qs (dirname, c, NESS = False, relax = False, saved_ts = None, saved_ps = None):
    '''
    Wrapper function used to generate and write any of the possible types of trajectories, probability distributions, etc.
    inputs already described in the functions that require them
    '''
    try:
        os.chdir(dirname)
    except ValueError:
        print ('Not a valid directory of saved trajectories!')

    '''
    Okay so after writing, now read the stored trajectory data and calculate quantities
    What's the best way to do this? It may perhaps be that these calculations are best done functionally instead of with a class -
    that is, if you have everything you need saved, then you can just load what you need to calculate; and calculate as you go

    you'll need to save those calculated quantities also as pytables files within the same directory as the trajectories
    NESS quantities
    Relaxation quantities

    To read the npy files:
    if dict of arrs (e.g., X_i):
        varname = np.load(filepath, allow_pickle=True)[()]
    else (e.g., X):
        varname = np.load(filepath, allow_pickle=True)
    '''

    coordinates = np.load('coordinates.npy', allow_pickle=True)
    num_coords = len(coordinates)
    mechanisms = np.load('mechanisms.npy', allow_pickle=True)
    num_mechs = len(mechanisms)

    puppets = np.load('puppets.npy', allow_pickle=True)[()]
    leaders = np.load('leaders.npy', allow_pickle=True)[()]
    negleaders = np.load('negleaders.npy', allow_pickle=True)[()]
    units = np.load('units.npy', allow_pickle=True)[()]
    num_units = len(units)

    X = np.load('X.npy', allow_pickle=True)
    X_i = np.load('X_i.npy', allow_pickle=True)[()]
    X_u = np.load('X_u.npy', allow_pickle=True)[()]
    X_p = np.load('X_p.npy', allow_pickle=True)[()]
    X_l = np.load('X_l.npy', allow_pickle=True)[()]
    X_nl = np.load('X_nl.npy', allow_pickle=True)[()]

    p_ness = np.load('p_ness.npy', allow_pickle=True)
    p_init = np.load('p_init.npy', allow_pickle=True)

    Ksys = np.load('Ksys.npy', allow_pickle=True)
    Ks_v = np.load('Ks_v.npy', allow_pickle=True)[()]

    sys_props = System_Properties(  coordinates,
                                    num_coords,
                                    mechanisms,
                                    num_mechs,
                                    puppets,
                                    leaders,
                                    negleaders,
                                    units,
                                    num_units,
                                    X,
                                    X_i,
                                    X_u,
                                    X_p,
                                    X_l,
                                    X_nl,
                                    p_ness,
                                    p_init,
                                    Ksys,
                                    Ks_v)

    class cols_qs(tb.IsDescription):
        time = tb.Float64Col(pos = 0)

        if_mechs = tb.Float64Col(pos = 1, shape = (num_mechs, ))
        if_coords = tb.Float64Col(pos = 2, shape = (2*choose_two(num_coords), ))
        if_puppets = tb.Float64Col(pos = 3, shape = (2*choose_two(num_mechs), ))
        if_leaders = tb.Float64Col(pos = 4, shape = (2*choose_two(num_mechs), ))
        if_units = tb.Float64Col(pos = 5, shape = (2*choose_two(num_units), ))

        mi_inst_mechs = tb.Float64Col(pos = 6, shape = (num_mechs, ))
        mi_inst_coords = tb.Float64Col(pos = 7, shape = (choose_two(num_coords), ))
        mi_inst_puppets = tb.Float64Col(pos = 8, shape = (choose_two(num_mechs), ))
        mi_inst_leaders = tb.Float64Col(pos = 9, shape = (choose_two(num_mechs), ))
        mi_inst_units = tb.Float64Col(pos = 10, shape = (choose_two(num_units), ))

        mi_tot_mechs = tb.Float64Col(pos = 11, shape = (num_mechs, ))
        mi_tot_coords = tb.Float64Col(pos = 12, shape = (choose_two(num_coords), ))
        mi_tot_puppets = tb.Float64Col(pos = 13, shape = (choose_two(num_mechs), ))
        mi_tot_leaders = tb.Float64Col(pos = 14, shape = (choose_two(num_mechs), ))
        mi_tot_units = tb.Float64Col(pos = 15, shape = (choose_two(num_units), ))

        ep = tb.Float64Col(pos = 16)
        zeta_mechs = tb.Float64Col(pos = 17, shape = (num_mechs, ))
        ep_units = tb.Float64Col(pos = 18, shape = (num_units, ))

        ef = tb.Float64Col(pos = 18)
        ef_mechs = tb.Float64Col(pos = 20, shape = (num_mechs, ))
        ef_puppets = tb.Float64Col(pos = 21, shape = (num_mechs, ))
        ef_leaders = tb.Float64Col(pos = 22, shape = (num_mechs, ))
        ef_units = tb.Float64Col(pos = 23, shape = (num_units, ))

        s = tb.Float64Col(pos = 24)
        s_coords = tb.Float64Col(pos = 25, shape = (num_coords, ))
        s_puppets = tb.Float64Col(pos = 26, shape = (num_mechs, ))
        s_leaders = tb.Float64Col(pos = 27, shape = (num_mechs, ))
        s_units = tb.Float64Col(pos = 28, shape = (num_units, ))


    if NESS:
        print('Writing NESS quantities...')
        calc_qs(c, '_traj_NESS.h5', '_traj_quants_NESS.h5', sys_props, cols_qs)
    if relax:
        print('Writing relaxation quantities...')
        calc_qs(c, '_traj_rlx.h5', '_traj_quants_rlx.h5', sys_props, cols_qs, rlx = True)
        # use saved_ts and saved_ps here

    os.chdir('..')


this_model_name = 'test_ckt'
this_dirname = this_model_name

# write_qs(this_dirname, this_model_name, NESS = True)
write_qs(this_dirname, this_model_name, relax = True)
