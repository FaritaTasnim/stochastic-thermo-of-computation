import numpy as np
import tables as tb
from itertools import product
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
rc('animation', html='html5')
from IPython.display import HTML
from scipy.sparse.linalg import eigs
from scipy.stats import pearsonr
import matplotlib as mpl
from matplotlib.colors import LogNorm, PowerNorm, Normalize
from matplotlib.lines import Line2D


font = {'weight': 'bold',
        'size'   : 16}

plt.rc('font', **font)


rng = np.random.default_rng()





'''
################################################################################
Setting all model conditions from which you want to generate plots
################################################################################
'''

n_traj = 200
n_wind = 200
n_steps = 2001

initdists = ['Cness']
energy_fcns = ['delta', 'abs']
# red red red
# remove below line for when re-collecting files
# energy_fcns = [r'$-\delta(x_i - x_o)$', r'$|x_i - x_o|$']
n_msgs = np.linspace(2, 5, 4).astype(int)
To_vals = np.round(np.linspace(0.5, 5.5, 11), 2)
R1_sr_vals = [0.0]
NR_sr_vals = np.round(np.linspace(0.5, 5.5, 11), 2)
steps = np.round(np.linspace(1, n_steps, n_steps), 2)


'''
################################################################################
Generating reference arrays with model conditions
################################################################################
'''

'''
Everything should be stored in an array with dims starting with
len(initdists), len(energy_fcns), len(n_msgs), len(To_vals), len(sr_vals)
'''

def create_refs(ids, efs, nms, Tovs, srvs):
    ref = np.empty((len(ids), len(efs), len(nms), len(Tovs), len(srvs), 5), dtype=object)
    for na, a in enumerate(ids):
        for nb, b in enumerate(efs):
            for nc, c in enumerate(nms):
                for nd, d in enumerate(Tovs):
                    for ne, e in enumerate(srvs):
                        ref[na, nb, nc, nd, ne] = [a,b,c,d,e]
    return ref

ref_R1 = create_refs(initdists, energy_fcns, n_msgs, To_vals, R1_sr_vals)
ref_NR = create_refs(initdists, energy_fcns, n_msgs, To_vals, NR_sr_vals)
refs = {'r1': ref_R1,
        'nr': ref_NR}

conds_dict = {  'r1': {'id': initdists, 'H': energy_fcns, 'L': n_msgs, 'To': To_vals, 'sr': R1_sr_vals},
                'nr': {'id': initdists, 'H': energy_fcns, 'L': n_msgs, 'To': To_vals, 'sr': NR_sr_vals}}



dims = {m: tuple([len(conds_dict[m][cond]) for cond in ['id', 'H', 'L', 'To', 'sr']]) for m in ['r1', 'nr']}

def extended_dims(start_dims, dims_to_add):
    end_dims = list(start_dims)
    for dim in dims_to_add:
        end_dims.append(dim)
    return tuple(end_dims)

dims_traj = {m: extended_dims(dims[m], [n_traj]) for m in ['r1', 'nr']}
dims_wind = {m: extended_dims(dims[m], [n_wind]) for m in ['r1', 'nr']}
dims_traj_rlx = {m: extended_dims(dims[m], [n_steps, n_traj]) for m in ['r1', 'nr']}
dims_fcn_rlx = {m: extended_dims(dims[m], [n_steps]) for m in ['r1', 'nr']}

dim_iters = {m: tuple([np.arange(d) for d in dims[m]]) for m in ['r1', 'nr']}
dims_fcn_rlx_iters = {m: tuple([np.arange(d) for d in dims_fcn_rlx[m]]) for m in ['r1', 'nr']}


'''
For example to get the list of conditions for everything specified but To_val
ref[0,0,1,:,0]
or for everything specified but To_val and sr_val
ref[0,0,1,:,:]
'''
# print(ref_R1[0,0,1,:,:])
# print()
# print(ref_NR[0,0,1,:,:])

base_axes = { 'id': 0, 'H': 1, 'L': 2, 'To': 3, 'sr': 4, 'num_axes': 5}
traj_axes = { 'id': 0, 'H': 1, 'L': 2, 'To': 3, 'sr': 4, 'traj': 5, 'num_axes': 6}
time_axes = { 'id': 0, 'H': 1, 'L': 2, 'To': 3, 'sr': 4, 'time': 5, 'num_axes': 6}
timextraj_axes = { 'id': 0, 'H': 1, 'L': 2, 'To': 3, 'sr': 4, 'time': 5, 'traj': 6, 'num_axes': 7}

qcolors = { 'if_io':    'fuchsia',
            'if_oi':    'orchid',
            'mi_inst':  'darkmagenta',
            'mi_tot':   'indigo',
            'ep':       'crimson',
            'q':        'orange',
            's':        'turquoise',
            'x':        'deepskyblue'}




'''
################################################################################
Getting filenames
################################################################################
'''

def get_filename(modelname, filetype, id, energy_fcn, num_msg, To_val, sr_val):
    this_model_name = 'model_{}_{}-{}_{}-{}_{}-{}_{}-{}_{}-{}'.format(modelname, 'H', energy_fcn, 'To', To_val, 'sr', sr_val, 'id', id, 'L', num_msg)
    return '{}_{}.h5'.format(this_model_name, filetype)

def get_all_filenames(modelname, p_or_t, n_or_r, t_or_w):
    '''
    p_or_t: Probability distributions or Trajectories,
    n_or_r: NESS or Relaxation,
    t_or_w: Trajectories or Windows
    '''

    fnames = np.empty(dims[modelname], dtype=object)

    if p_or_t == 'p':
        filetype = 'pdists_relax'
    else:
        if n_or_r == 'r':
            filetype = 'quants_relax'
        else:
            if t_or_w == 't':
                filetype = 'traj_quants_NESS'
            elif t_or_w == 'w':
                filetype = 'wind_quants_NESS'

    for na, nb, nc, nd, ne in product(*dim_iters[modelname]):
        vals = refs[modelname][na,nb,nc,nd,ne]
        fnames[na,nb,nc,nd,ne] = get_filename(modelname, filetype, *vals)

    return fnames

r1_pdist_rlx_qfs = get_all_filenames('r1', 'p', '-', '-')
r1_traj_rlx_qfs = get_all_filenames('r1', '-', 'r', '-')
r1_traj_ness_qfs = get_all_filenames('r1', '-', 'n', 't')

nr_pdist_rlx_qfs = get_all_filenames('nr', 'p', '-', '-')
nr_traj_rlx_qfs = get_all_filenames('nr', '-', 'r', '-')
nr_traj_ness_qfs = get_all_filenames('nr', '-', 'n', 't')
nr_wind_ness_qfs = get_all_filenames('nr', '-', 'n', 'w')
print('Got all filenames!')
print()

# print(nr_wind_ness_qfs[0,0,:,:,:])







'''
################################################################################
Calculation helper functions
################################################################################
'''

'''Note: Have to change this for a chain'''
input_sub = 0
output_sub = 1

def get_p_unequal(p_io, attrs, axis):
    ss_io = attrs.space
    where_unequal = np.where(ss_io[:,input_sub] != ss_io[:,output_sub])[0]
    p_unequals = np.take(p_io, where_unequal, axis = axis)
    # print(len(p_unequals.shape))
    if len(p_unequals.shape) == 1:
        return np.sum(p_unequals, axis=0)
    else:
        return np.sum(p_unequals, axis=1)







'''
################################################################################
Getting quantities or pdists from files
################################################################################
'''

np.set_printoptions(suppress=True, precision=3)

def get_all_Cness (modelname, fs):
    all_Cness = np.empty(fs.shape)
    for i in np.ndindex(fs.shape):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_Cness[i] = f.root._v_attrs.Cness
        f.close()
    return all_Cness

try:
    r1_Cness = np.load('r1_Cness.npy')
except:
    print('No R1 Cness values saved. Regenerating...')
    r1_Cness = get_all_Cness('r1', r1_traj_ness_qfs)
    np.save('r1_Cness.npy', r1_Cness)
    print('Got em!')

try:
    nr_Cness = np.load('nr_Cness.npy')
except:
    print('No NR Cness values saved. Regenerating...')
    nr_Cness = get_all_Cness('nr', nr_traj_ness_qfs)
    np.save('nr_Cness.npy', nr_Cness)
    print('Got em!')

# r1_Cness = get_all_Cness('r1', r1_traj_ness_qfs)
# nr_Cness = get_all_Cness('nr', nr_traj_ness_qfs)

Cness = {'r1': r1_Cness, 'nr': nr_Cness}

print('Got all Cness values!')
print()

# print(nr_Cness)
# print(refs['nr'][0,0,:,:,:,2])
# print(refs['nr'][0,0,:,:,:,3])
# print(refs['nr'][0,0,:,:,:,4])

def get_all_rlx_times (modelname, fs):
    all_rlx_times = np.empty(dims_fcn_rlx[modelname]) # dims x num_steps
    for i in np.ndindex(fs.shape):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_rlx_times[i] = f.root.pdists_relax.pdists_relax.col('time')[:n_steps]
        f.close()
    return all_rlx_times

# r1_rlx_times = get_all_rlx_times('r1', r1_pdist_rlx_qfs)
# nr_rlx_times = get_all_rlx_times('nr', nr_pdist_rlx_qfs)

try:
    r1_rlx_times = np.load('r1_rlx_times.npy')
except:
    print('No R1 relaxation times saved. Regenerating...')
    r1_rlx_times = get_all_rlx_times('r1', r1_pdist_rlx_qfs)
    np.save('r1_rlx_times.npy', r1_rlx_times)
    print('Got em!')

try:
    nr_rlx_times = np.load('nr_rlx_times.npy')
except:
    print('No NR relaxation times saved. Regenerating...')
    nr_rlx_times = get_all_rlx_times('nr', nr_pdist_rlx_qfs)
    np.save('nr_rlx_times.npy', nr_rlx_times)
    print('Got em!')

times = {'r1': r1_rlx_times, 'nr': nr_rlx_times}

print('Got all relaxation times!')
print()

# print(r1_rlx_times)
# print(r1_rlx_times.shape)

def get_all_ness_quants (modelname, fs, quant, dims_t=dims_traj):
    ''' Can change dims_traj to dims_wind if you decide to increase number of windows collected'''
    all_quants = np.empty(dims_t[modelname]) # dims x num_traj or num_wind
    for i in tqdm(np.ndindex(fs.shape)):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_quants[i] = f.root.quants.quants.col(quant)
        f.close()
    return all_quants

# r1_ep_traj_ness = get_all_ness_quants('r1', r1_traj_ness_qfs, 'ep')
# nr_ep_wind_ness = get_all_ness_quants('nr', nr_wind_ness_qfs, 'ep')

# print(nr_ep_wind_ness[0,0,:,:,:])
# print(nr_ep_wind_ness.shape)


def get_all_rlx_quants (modelname, fs, quant):
    all_rlx_quants = np.empty(dims_traj_rlx[modelname]) # dims x num_steps x num_traj
    for i in tqdm(np.ndindex(fs.shape)):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_rlx_quants[i] = np.concatenate([[run.col(quant)] for run in f.root.quants], axis = 0).T[:n_steps,:]
        f.close()
    return all_rlx_quants

# nr_ifio_traj_rlx = get_all_rlx_quants('nr', nr_traj_rlx_qfs, 'if_io')
# print(nr_ifio_traj_rlx[0,0,:,:,:])
# print(nr_ifio_traj_rlx.shape)

def get_all_fin_rlx_quants (modelname, fs, quant):
    all_fin_rlx_quants = np.empty(dims_traj[modelname]) # dims x num_traj
    for i in tqdm(np.ndindex(fs.shape)):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_fin_rlx_quants[i] = [run.col(quant)[-1] for run in f.root.quants]
        f.close()
    return all_fin_rlx_quants

# nr_ifio_fin_rlx = get_all_fin_rlx_quants('nr', nr_traj_rlx_qfs, 'if_io')
# print(nr_ifio_fin_rlx[0,0,:,:,:])
# print(nr_ifio_fin_rlx.shape)

def get_all_ness_fcn_pios (modelname, fs, fcn):
    all_fcn_pios = np.empty(dims[modelname])
    for i in tqdm(np.ndindex(fs.shape)):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_fcn_pios[i] = fcn(f.root._v_attrs.p_io, f.root._v_attrs, axis=0)
        f.close()
    return all_fcn_pios

# nr_pness_unequal = get_all_ness_fcn_pios('nr', nr_traj_ness_qfs, get_p_unequal)
# print(nr_pness_unequal[:,:,:,:,:])
# print(nr_pness_unequal.shape)

def get_all_rlx_fcn_pios (modelname, fs, fcn):
    all_rlx_fcn_pios = np.empty(dims_fcn_rlx[modelname]) # dims x num_steps
    for i in tqdm(np.ndindex(fs.shape)):
        fname = fs[i]
        f = tb.open_file(fname, mode='r')
        all_rlx_fcn_pios[i] = fcn(  f.root.pdists_relax.pdists_relax.col('p_io')[:,0],
                                    f.root.pdists_relax._v_attrs, axis=1)[:n_steps]
        f.close()
    return all_rlx_fcn_pios

# nr_prelax_unequal = get_all_rlx_fcn_pios('nr', nr_pdist_rlx_qfs, get_p_unequal)
# print(nr_prelax_unequal[:,:,:,:,:])
# print(nr_prelax_unequal.shape)







'''
################################################################################
Functions of quantities
################################################################################
'''
def no_fcn (qs):
    return qs

def avg (qs):
    return np.mean(qs, axis = -1)

def avg_reciprocal (qs):
    return np.mean(1 / qs, axis = -1)

def avg_twice_reciprocal (qs):
    return np.mean(2 / qs, axis = -1)

def reciprocal_avg (qs):
    return 1 / np.mean(qs, axis = -1)

def twice_reciprocal_avg (qs):
    return 2*reciprocal_avg(qs)

def lack_of_precision (qs):
    return np.var(qs, axis = -1) / (np.mean(qs, axis = -1)**2)

def precision (qs):
    return (np.mean(qs, axis = -1)**2) / np.var(qs, axis = -1)

def variance (qs):
    return np.var(qs, axis = -1)

def ift (eps):
    return np.mean(np.exp(-1*eps), axis = -1)

def remove_nan (a):
    a[np.isnan(a)] = 0.

def cov (qpairs):
    all_covs = np.zeros(qpairs.shape[:-2])
    for i in np.ndindex(qpairs.shape[:-2]):
        all_covs[i] = np.corrcoef(qpairs[i], rowvar=False)[0,1]**2
    remove_nan(all_covs)
    return all_covs

def pearson_cc (qpairs):
    all_pccs = np.zeros(qpairs.shape[:-2])
    for i in np.ndindex(qpairs.shape[:-2]):
        r, _ = pearsonr(qpairs[i][...,0] , qpairs[i][...,1])
        all_pccs[i] = r
    remove_nan(all_pccs)
    return all_pccs

def squared_pearson_cc (qpairs):
    return (pearson_cc(qpairs))**2


# def cov_time(modelname, qpairs):
#     all_covs = np.empty(dims_fcn_rlx[modelname]) # dims x num_steps
#     for na, nb, nc, nd, ne, nf in product(*dims_fcn_rlx_iters[modelname]):
#         all_covs[na,nb,nc,nd,ne,nf] = np.corrcoef(qpairs[na,nb,nc,nd,ne,nf,:], rowvar=False)[0,1]**2
#     remove_nan(all_covs)
#     return all_covs


'''
For functions of multiple quantities
qpairs shape:
    dims x num_steps x num_traj x num_quants, or
    dims x num_traj x num_quants
'''

# nr_ifio_traj_rlx = get_all_rlx_quants('nr', nr_traj_rlx_qfs, 'if_io')
# nr_mitot_traj_rlx = get_all_rlx_quants('nr', nr_traj_rlx_qfs, 'mi_tot')
# nr_ifiomitotpairs_traj_rlx = np.concatenate([np.expand_dims(nr_ifio_traj_rlx, -1), np.expand_dims(nr_mitot_traj_rlx, -1)], axis=-1)
# covs = cov_time('nr', nr_ifiomitotpairs_traj_rlx)
#
# print(nr_ifiomitotpairs_traj_rlx[0,0,0,0,0,:].shape)
# print(covs)
# print(covs.shape)





'''
################################################################################
Plotting helper functions
################################################################################
'''

def downsample_along_axis (arr, num, axis, start=0):
    # print(arr.shape)
    indices = np.linspace(start, arr.shape[axis] - 1, num).astype(int)
    return np.take(arr, indices, axis = axis)

downsampled_t = lambda t, dsamp: downsample_along_axis(t, dsamp, -1) if dsamp is not None else t

def get_quants (q_reaper, m, qfs, q, qfname, dsample, axes, which_axis):
    qf_npy = '{}_{}.npy'.format(qfname, q)
    print(qf_npy)
    try:
        vals = np.load(qf_npy)
        print('Loaded saved quants file!')
    except:
        print('{} not saved. Generating...'.format(qf_npy))
        vals = q_reaper(m, qfs, q)
        np.save(qf_npy, vals)
        print('Got em!')
    # vals = q_reaper(m, qfs, q)
    if dsample is not None:
        if q == 'ep' or q == 's':
            vals = downsample_along_axis(vals, dsample, axes[which_axis], start=1)
        else:
            vals = downsample_along_axis(vals, dsample, axes[which_axis])
        remove_nan(vals)
    return vals

def apply_fcn_to_quants (   q_reaper, m, qfs, qs, qfname, fcn,
                            dsample, axes, x_axis):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    qs (str) or (list): name(s) of quantity to extract,
    fcn (fcn): to apply to quantit(ies)
    dsample (float): if provided, how much to downsample x_axis,
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): along which axis to make line graph
    '''
    if isinstance(qs, list):
        all_qs = [get_quants(q_reaper, m, qfs, q, qfname, dsample, axes, x_axis) for q in qs]
        vals = np.concatenate([np.expand_dims(all_q, -1) for all_q in all_qs], axis=-1)
    else:
        vals = get_quants(q_reaper, m, qfs, qs, qfname, dsample, axes, x_axis)
    return fcn(vals)

def get_vals_conds( q_reaper, m, qfs, qs, qfname, fcn,
                    axes, x_axis, dsample, z_axis):

    fd_vals = apply_fcn_to_quants(q_reaper, m, qfs, qs, qfname, fcn, dsample, axes, x_axis)
    if x_axis != 'time':
        fd_vals, new_axes = roll_vals_axes(fd_vals, axes, x_axis, axes['num_axes'])
        fd_vals, new_axes = roll_vals_axes(fd_vals, new_axes, z_axis, axes['num_axes']-1)
        new_times = None
    else:
        new_times, _ = roll_vals_axes(times[m], axes, z_axis, axes['num_axes']-1)
        fd_vals, new_axes = roll_vals_axes(fd_vals, axes, z_axis, axes['num_axes']-1)

    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[x_axis], new_axes[z_axis], new_axes['num_axes']])
    if x_axis != 'time':
        xvals = conds_dict[m][x_axis]
    else:
        xvals = None
    zvals = conds_dict[m][z_axis]

    return fd_vals, conds, xvals, zvals, new_times

roll_axes = lambda axes, x_axis, dest: {key: val if key == 'num_axes' else val-1 if val > axes[x_axis] and val <= dest-1 else dest-1 if key == x_axis else val for i, (key, val) in enumerate(axes.items())}

def roll_vals_axes (vals, axes, which_axis, destination):
    vals = np.rollaxis(vals, axes[which_axis], destination)
    new_axes = roll_axes(axes, which_axis, destination)
    return vals, new_axes


def plot_2d_lines ( toplot, xvals, zvals,
                    xlab, ylab, zlab, zcolors,
                    xlog, ylog, savename, ymax=None, ymin=None):

    plt.figure(figsize=(8, 8), dpi=300)

    for i, (data, zval) in enumerate(zip(toplot, zvals)):
        if len(xvals.shape) > 1:
            plt.plot(xvals[i], data, linewidth = 7, alpha = 0.8, color = zcolors[i], label = '{} = {}'.format(zlab, zval))
        else:
            plt.plot(xvals, data, linewidth = 7, alpha = 0.8, color = zcolors[i], label = '{} = {}'.format(zlab, zval))

    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    if ymax is not None:
        plt.ylim(top=ymax)
    if ymin is not None:
        plt.ylim(bottom=ymin)

    if ylog:
        plt.yscale('log')
    if xlog:
        plt.xscale('log')

    plt.tight_layout()
    plt.savefig(savename)
    plt.close()


def discrete_cmap(base_cmap, N):
    base = mpl.cm.get_cmap(base_cmap)
    if N > 1:
        color_list = np.linspace(0, 1, N, 0)
        color_list = color_list + color_list[1]
        return color_list, base, [base(c) for c in color_list]
    else:
        color_list = [1.]
        return color_list, base, [base(c) for c in color_list]

markers = ['o', 'v', 's', '*', '+', 'p', '^', 'x', '1', 'd', 'D']
qmarkers = ['d', 'P', '^', 'o', 's']
pcolors = { 'id': discrete_cmap('Greys', len(initdists)),
            'H': discrete_cmap('YlOrBr', len(energy_fcns)),
            'L': discrete_cmap('RdPu', len(n_msgs)),
            'To': discrete_cmap('coolwarm', len(To_vals)),
            'sr_r1': discrete_cmap('plasma_r', len(R1_sr_vals)),
            'sr_nr': discrete_cmap('viridis', len(NR_sr_vals))}

def plot_2d_scatter(qs, xs, qlab, xlab, cs, ms,
                    clabs, mlabs, savename, ylog,
                    mult=False, cname='', mname='q',
                    ymax=None, ymin=None, legendloc=None):

    plt.figure(figsize=(8, 8), dpi=300)
    ax = plt.gca()
    cinds, cmap, _ = cs
    if mult:
        Cind, M = np.meshgrid(cinds, ms, indexing='xy')
    else:
        Cind, M = np.meshgrid(cinds, ms, indexing='ij')

    if len(Cind.shape) == len(xs.shape) - 1:
        Cind = np.expand_dims(Cind, -1)
        Cind = np.repeat(Cind, xs.shape[-1], -1)
        M = np.expand_dims(M, -1)
        M = np.repeat(M, xs.shape[-1], -1)

    for i in np.ndindex(xs.shape):
        ax.scatter(xs[i], qs[i], color=cmap(Cind[i]), marker=M[i], s=100)

    if ylog:
        ax.set_yscale('log')

    clegend = [Line2D([0], [0], marker='o', color=cmap(cinds[i]), lw=0, label='{} = {}'.format(cname,clab), markersize=8) for i, clab in enumerate(clabs)]
    mlegend = [Line2D([0], [0], marker=ms[i], color='gray', lw=0, label='{} = {}'.format(mname,mlab), markersize=8) for i, mlab in enumerate(mlabs)]
    legend_elems = clegend + mlegend

    if legendloc is not None:
        plt.legend(handles=legend_elems, loc=legendloc)
    else:
        plt.legend(handles=legend_elems)

    plt.xlabel(xlab)
    plt.ylabel(qlab)

    if ymax is not None:
        plt.ylim(top=ymax)
    if ymin is not None:
        plt.ylim(bottom=ymin)


    plt.tight_layout()
    plt.savefig(savename)
    plt.close()

def animate_histogram (all_vals, num_bins, col, fps, animsavename, xlab):
    # shape of all_vals: num_frames (which corresponds either to parameter values or time values) x num_traj
    qmin, qmax = np.amin(all_vals), np.amax(all_vals)
    bin_list = np.linspace(qmin, qmax, num_bins + 1)

    def prepare_animation(bar_container):
        def animate(frame_number):
            data = all_vals[frame_number]
            n, _ = np.histogram(data, bin_list)
            for count, rect in zip(n, bar_container):
                rect.set_height(count)
            return bar_container
        return animate

    fig, ax = plt.subplots()
    _, _, bar_container = ax.hist(all_vals[0], bin_list, lw=1, ec='white', fc=col, alpha=0.5)
    ax.set_ylim(top=all_vals.shape[1]/2)  # set safe limit to ensure that all data is visible.
    ax.set_xlabel(xlab)

    ani = animation.FuncAnimation(fig, prepare_animation(bar_container), all_vals.shape[0], repeat=False, blit=True)

    # plt.show()
    writergif = animation.PillowWriter(fps=fps)
    ani.save(animsavename, writer=writergif)
    plt.close()


def heatmap (data, xlab, ylab, row_labels=None, col_labels=None, title=None, ax=None, cbar_kw={}, cbarlabel="", logscale=False, powerscale=False, **kwargs):

    # Plot the heatmap
    if logscale:
        im = ax.imshow(data, norm=LogNorm(), **kwargs)
    elif powerscale:
        im = ax.imshow(data, norm=PowerNorm(0.3), **kwargs)
    else:
        im = ax.imshow(data, **kwargs)

    # Create colorbar
    if logscale or powerscale:
        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.75, format='%.e', **cbar_kw)
    else:
        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.75, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    num_xticks = int(data.shape[1]/5) + 1
    num_yticks = int(data.shape[0]/3) + 1
    # print(num_xticks, num_yticks)
    ax.set_xticks(np.arange(num_xticks)*5)
    ax.set_yticks(np.arange(num_yticks)*3)
    if col_labels is not None:
        print(col_labels.shape)
        print(np.arange(num_xticks)*5)
        ax.set_xticklabels(col_labels[np.arange(num_xticks)*5])
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
    if row_labels is not None:
        ax.set_yticklabels(row_labels[np.arange(num_yticks)*3])

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Turn spines off and create white grid.
    for key, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    if title != None:
        ax.set_title(title)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar




'''
################################################################################
S1, Cness vs params, 2D
################################################################################
'''

def plot_Cness_vparams( m, x_axis, m_axis, c_axis,
                        ms, cs, xname, mname, cname, t_or_w,
                        Cnesslog, legloc):
    '''
    m (str): model name,
    m_axis (str): one marker style for each of this parameter's values,
    c_axis (str): one color style for each of this parameter's values,
    ms (list) one for each value m_axis takes on,
    cs (list) one for each value c_axis takes on,
    qlab (str): label for quantity being plotted on y-axis,
    mlab (str): label for m_axis,
    clab (str): label for c_axis,
    t_or_w (str): 't' if for trajectories, 'w' if for windows
    '''

    # print(Cness[m].shape)
    new_Cness, new_axes = roll_vals_axes(Cness[m], base_axes, x_axis, base_axes['num_axes'])
    xvals, na = roll_vals_axes(refs[m][...,base_axes[x_axis]], base_axes, x_axis, base_axes['num_axes'])
    new_Cness, new_axes = roll_vals_axes(new_Cness, na, m_axis, base_axes['num_axes']-1)
    xvals, na = roll_vals_axes(xvals, na, m_axis, base_axes['num_axes']-1)
    new_Cness, new_axes = roll_vals_axes(new_Cness, new_axes, c_axis, base_axes['num_axes']-2)
    xvals, na = roll_vals_axes(xvals, na, c_axis, base_axes['num_axes']-2)

    # can function-ify this red red red
    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[x_axis], new_axes[m_axis], new_axes[c_axis], new_axes['num_axes']])

    marks = ms[:len(conds_dict[m][m_axis])]
    clabs = conds_dict[m][c_axis]
    mlabs = conds_dict[m][m_axis]


    for i in np.ndindex(new_Cness.shape[:-3]):
        savename = 'model_{}_{}_Cness_vs-{}_for-{}s-and-{}s_{}.png'.format(m, t_or_w, x_axis, c_axis, m_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        plot_2d_scatter(new_Cness[i], xvals[i], r'$C^{NESS}$', xname, cs, marks,
                        clabs, mlabs, savename, Cnesslog, mult=0, cname=cname, mname=mname,
                        ymax=1.1*new_Cness.max(), ymin=0.9*new_Cness.min(), legendloc= legloc)



'''
################################################################################
A1, Animated histograms, 2D
################################################################################
'''

'''A1, Relaxation'''
'''A1, NESS'''
def animate_2d_dists(   q_reaper, m, qfs, q, qfname,
                        axes, anim_axis, dsample,
                        nbins, c, fps, qlab, t_or_w,
                        process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    q (str): name of quantity to extract,
    axes (dict): of 'this': axis along which 'this' varies,
    anim_axis (str): along which axis to animate,
    dsample (float): if provided, how much to downsample anim_axis,
    nbins (int): number of bins for histogram,
    c (str): color for histogram bars,
    fps (int): frames per second for animation,
    qlab (str): label for x-axis of histogram,
    t_or_w (str): 't' if for trajectories, 'w' if for windows
    process (str): identifying the process, e.g. 'rlx' or 'fin_rlx' or 'ness'
    '''
    vals = get_quants(q_reaper, m, qfs, q, qfname, dsample, axes, anim_axis)
    vals, new_axes = roll_vals_axes(vals, axes, anim_axis, axes['num_axes']-1)
    # np.rollaxis(vals, axes[anim_axis], -1) # last two axes: anim_axis x trajs

    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[anim_axis], new_axes['traj'], new_axes['num_axes']])

    for i in np.ndindex(vals.shape[:-2]):
        savename = 'model_{}_{}_{}_anim-{}_vs-{}_{}.gif'.format(m, process, t_or_w, q, anim_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        animate_histogram(vals[i], nbins, c, fps, savename, qlab)






'''
################################################################################
S2, Static plot of f(quants) vs param, 2D
################################################################################
'''

def plot_2d_fcn_of_quants(  q_reaper, m, qfs, qs, fcn,
                            axes, x_axis, dsample, z_axis,
                            xlab, qlab, zlab, zcolors, t_or_w,
                            xlog, ylog, savelab, process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    qs (str) or (list): name(s) of quantity to extract,
    fcn (fcn): to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): along which axis to make line graph,
    dsample (float): if provided, how much to downsample x_axis,
    z_axis (str): one line for each of this parameter's values,
    xlab (str): label for parameter on x-axis,
    qlab (str): label for quantity being plotted on y-axis,
    zlab (str): label for parameter associated with each line,
    zcolors (list): of strings giving colors for line pertaining to each z,
    t_or_w (str): 't' if for trajectories, 'w' if for windows,
    xlog (bool): True if x-axis should be logscaled,
    ylog (bool): True if y-axis should be logscaled,
    savelab (str): identifier of this plot in the name of the png,
    '''
    fd_vals, conds, xvals, zvals, new_times = get_vals_conds(   q_reaper, m, qfs, qs, fcn,
                                                                axes, x_axis, dsample, z_axis)
    _, _, colors = pcolors[z_axis]
    for i in np.ndindex(fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-{}_for-{}s_{}.png'.format(m, process, t_or_w, savelab, x_axis, z_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        if x_axis == 'time':
            xvals = downsampled_t(new_times, dsample)[i]
        plot_2d_lines(fd_vals[i], xvals, zvals, xlab, qlab, zlab, colors, xlog, ylog, savename)


# plt.rc('legend',fontsize='small')

# plot_2d_fcn_of_quants(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_tot'], cov,
#                         base_axes, x_axis='To', dsample=None, z_axis='L',
#                         xlab=r'$T^o_2$',
#                         qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                         zlab='L', zcolors='', t_or_w='w',
#                         xlog=0, ylog=0, savelab='ifio-mitot-corr')
# #
# plot_2d_fcn_of_quants(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'], cov,
#                         time_axes, x_axis='time', dsample=51, z_axis='L',
#                         xlab=r'$t$',
#                         qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                         zlab='L', zcolors='', t_or_w='t',
#                         xlog=0, ylog=0, savelab='ifio-mitot-corr')


def plot_2d_fcn_of_quants_vCness(   q_reaper, m, qfs, qs, qfname, fcn,
                                    axes, m_axis, c_axis,
                                    ms, cs, qlab, mname, cname, t_or_w,
                                    ylog, savelab, process='', legloc=None):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    qs (str) or (list): name(s) of quantity to extract,
    fcn (fcn): to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    m_axis (str): one marker style for each of this parameter's values,
    c_axis (str): one color style for each of this parameter's values,
    ms (list) one for each value m_axis takes on,
    cs (list) one for each value c_axis takes on,
    qlab (str): label for quantity being plotted on y-axis,
    mlab (str): label for m_axis,
    clab (str): label for c_axis,
    t_or_w (str): 't' if for trajectories, 'w' if for windows
    '''

    fd_vals = apply_fcn_to_quants(q_reaper, m, qfs, qs, qfname, fcn, None, None, None)

    fd_vals, new_axes = roll_vals_axes(fd_vals, axes, m_axis, axes['num_axes'])
    new_Cness, na = roll_vals_axes(Cness[m], axes, m_axis, axes['num_axes'])
    fd_vals, new_axes = roll_vals_axes(fd_vals, new_axes, c_axis, axes['num_axes']-1)
    new_Cness, na = roll_vals_axes(new_Cness, na, c_axis, axes['num_axes']-1)

    xlab = r'$C^{NESS}$'
    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[m_axis], new_axes[c_axis], new_axes['num_axes']])
    marks = ms[:len(conds_dict[m][m_axis])]
    clabs = conds_dict[m][c_axis]
    mlabs = conds_dict[m][m_axis]
    for i in np.ndindex(fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-Cness_for-{}s-and-{}s_{}.png'.format(m, process, t_or_w, savelab, c_axis, m_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        plot_2d_scatter(fd_vals[i], new_Cness[i], qlab, xlab, cs, marks,
                        clabs, mlabs, savename, ylog, cname=cname, mname=mname)
                        # ymax=1.1*fd_vals.max(), ymin=0.9*fd_vals.min(), legendloc= legloc)

# plot_2d_fcn_of_quants_vCness(get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 ['if_io', 'mi_tot'], cov,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr'],
#                                 qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='ifio-mitot-corr')

# plot_2d_fcn_of_quants_vCness(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 ['if_io', 'mi_tot'], cov,
#                                 base_axes, m_axis='H', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$H$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='ifio-mitot-corr')





'''
################################################################################
S3, Static plot of multiple f(quants) vs param, 2D
################################################################################
'''

def plot_2d_mult_fcns_of_quants(    q_reaper, m, qfs, all_qs, qfname, all_fcns,
                                    axes, x_axis, dsample,
                                    xlab, qlabs, colors, t_or_w,
                                    xlog, ylog, savelab, process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    all_qs (list): of strings and lists with name(s) of quantity to extract,
    all_fcns (list): of fcns to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): along which axis to make line graph,
    dsample (float): if provided, how much to downsample x_axis,
    xlab (str): label for parameter on x-axis,
    qlabs (list): of strings, each a label for quantity being plotted on y-axis,
    colors (list): of strings giving colors for line pertaining to each fcn(qs),
    t_or_w (str): 't' if for trajectories, 'w' if for windows,
    xlog (bool): True if x-axis should be logscaled,
    ylog (bool): True if y-axis should be logscaled,
    savelab (str): identifier of this plot in the name of the png,
    '''

    all_fd_vals = [apply_fcn_to_quants(q_reaper, m, qfs, qs, qfname, fcn, dsample, axes, x_axis) for qs, fcn in zip(all_qs, all_fcns)]
    all_fd_vals = np.concatenate([np.expand_dims(v, -1) for v in all_fd_vals], axis=-1)
    all_fd_vals, new_axes = roll_vals_axes(all_fd_vals, axes, x_axis, axes['num_axes']+1)

    if x_axis == 'time':
        new_times = np.concatenate([np.expand_dims(times[m], -1) for f in all_fcns], axis=-1)
        new_times, _ = roll_vals_axes(new_times, axes, x_axis, axes['num_axes']+1)

    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[x_axis], new_axes['num_axes']])

    if x_axis != 'time':
        xvals = conds_dict[m][x_axis]

    for i in np.ndindex(all_fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-{}_for-{}.png'.format(m, process, t_or_w, savelab, x_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        if x_axis == 'time':
            xvals = downsampled_t(new_times, dsample)[i]
        plot_2d_lines(all_fd_vals[i], xvals, qlabs, xlab, '', '', colors, xlog, ylog, savename)


# plot_2d_mult_fcns_of_quants(get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['ep', ['if_io', 'mi_tot']], [avg, cov],
#                             base_axes, x_axis='To', dsample=None,
#                             xlab=r'$T^o_2$', qlabs=[r'$<\sigma>$', r'$\chi(I^{i \to o}, I^{tot}(i;o))$'],
#                             colors=[qcolors['ep'], qcolors['x']], t_or_w='w', xlog=0, ylog=0, savelab='stuff')

# plot_2d_mult_fcns_of_quants(get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['ep', ['if_io', 'mi_tot']], [avg, cov],
#                             time_axes, x_axis='time', dsample=51,
#                             xlab=r'$t$', qlabs=[r'$<\sigma>$', r'$\chi(I^{i \to o}, I^{tot}(i;o))$'],
#                             colors=[qcolors['ep'], qcolors['x']], t_or_w='w', xlog=0, ylog=0, savelab='stuff')


def plot_2d_mult_fcn_of_quants_vCness(  q_reaper, m, qfs, all_qs, qfname, all_fcns,
                                        axes, c_axis, ms, cs,
                                        fqlabs, cname, t_or_w,
                                        ylog, savelab, process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    all_qs (list): of strings and lists with name(s) of quantity to extract,
    all_fcns (list): of fcns to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    c_axis (str): one color style for each of this parameter's values,
    ms (list) one for each value m_axis takes on,
    cs (list) one for each value c_axis takes on,
    qlab (str): label for quantity being plotted on y-axis,
    mlab (str): label for m_axis,
    clab (str): label for c_axis,
    t_or_w (str): 't' if for trajectories, 'w' if for windows
    '''
    all_fd_vals = [apply_fcn_to_quants(q_reaper, m, qfs, qs, qfname, fcn, None, None, None) for qs, fcn in zip(all_qs, all_fcns)]
    all_fd_vals = np.concatenate([np.expand_dims(v, -1) for v in all_fd_vals], axis=-1)
    all_fd_vals, new_axes = roll_vals_axes(all_fd_vals, axes, c_axis, axes['num_axes'])

    new_Cness = np.concatenate([np.expand_dims(Cness[m], -1) for f in all_fcns], axis=-1)
    new_Cness, _ = roll_vals_axes(new_Cness, axes, c_axis, axes['num_axes'])

    xlab = r'$C^{NESS}$'
    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[c_axis], new_axes['num_axes']])
    marks = ms[:len(all_fcns)]
    clabs = conds_dict[m][c_axis]
    for i in np.ndindex(all_fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-Cness_for-{}s_{}.png'.format(m, process, t_or_w, savelab, c_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        plot_2d_scatter(all_fd_vals[i], new_Cness[i], '', xlab, cs, marks,
                        clabs, fqlabs, savename, ylog, mult=0, cname=cname)

# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', ['if_io', 'mi_tot']], [avg, cov],
#                                     base_axes, c_axis='sr',
#                                     ms=markers, cs=pcolors['sr'],
#                                     fqlabs=['A','B'], cname=r'$f_s$', t_or_w = 'w',
#                                     ylog=0, savelab='stuff')






'''
################################################################################
S4, Static plot heatmapping f(quants) vs param and time, 2D
################################################################################
'''

def plot_fcn_of_quants_heatmap( q_reaper, m, qfs, qs, qfname, fcn,
                                axes, x_axis, dsample, z_axis,
                                xlab, qlab, zlab, cmap, t_or_w,
                                logscale, savelab, process='', powerscale=0):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    qs (str) or (list): name(s) of quantity to extract,
    fcn (fcn): to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): along which axis to make line graph,
    dsample (float): if provided, how much to downsample x_axis,
    z_axis (str): one line for each of this parameter's values,
    xlab (str): label for parameter on x-axis,
    qlab (str): label for quantity being plotted on y-axis,
    zlab (str): label for parameter associated with each line,
    zcolors (list): of strings giving colors for line pertaining to each z,
    t_or_w (str): 't' if for trajectories, 'w' if for windows,
    xlog (bool): True if x-axis should be logscaled,
    ylog (bool): True if y-axis should be logscaled,
    savelab (str): identifier of this plot in the name of the png,
    '''
    fd_vals, conds, xvals, zvals, new_times = get_vals_conds(   q_reaper, m, qfs, qs, qfname, fcn,
                                                                axes, x_axis, dsample, z_axis)
    for i in np.ndindex(fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-{}-and-{}-heatmap_{}.png'.format(m, process, t_or_w, savelab, x_axis, z_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)
        if x_axis == 'time':
            xvals = downsampled_t(new_times, dsample)[i]
        toplot = fd_vals[i]
        fig, ax = plt.subplots(figsize=(toplot.shape[1]/5 + 5, toplot.shape[0]/5 + 3), dpi=300)
        im, _ = heatmap(toplot, xlab, zlab, row_labels = zvals,
                        ax=ax, cmap=cmap, cbarlabel=qlab, logscale=logscale, powerscale=powerscale)
        # plt.tight_layout()
        plt.savefig(savename)
        plt.close()


fcmaps = {  'mean':         'viridis',
            'precision':    'plasma',
            'pearson_cc':   'Greens'}


# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_tot'], cov,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab='L', cmap='viridis', t_or_w='w',
#                             logscale=0, savelab='ifio-mitot-corr')

# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'], cov,
#                             time_axes, x_axis='time', dsample=51, z_axis='L',
#                             xlab=r'$t$',
#                             qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab='L', cmap='viridis', t_or_w='t',
#                             logscale=0, savelab='ifio-mitot-corr')








'''
################################################################################
S5, Static plot surface of f(quants) vs param and time, 3D
################################################################################
'''

def plot_fcn_of_quants_3D(  q_reaper, m, qfs, qs, qfname, fcn,
                            axes, x_axis, dsample, z_axis,
                            xlab, qlab, zlab, cmap, t_or_w,
                            logscale, savelab, process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    qs (str) or (list): name(s) of quantity to extract,
    fcn (fcn): to apply to quantit(ies)
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): along which axis to make line graph,
    dsample (float): if provided, how much to downsample x_axis,
    z_axis (str): one line for each of this parameter's values,
    xlab (str): label for parameter on x-axis,
    qlab (str): label for quantity being plotted on y-axis,
    zlab (str): label for parameter associated with each line,
    zcolors (list): of strings giving colors for line pertaining to each z,
    t_or_w (str): 't' if for trajectories, 'w' if for windows,
    xlog (bool): True if x-axis should be logscaled,
    ylog (bool): True if y-axis should be logscaled,
    savelab (str): identifier of this plot in the name of the png,
    '''
    fd_vals, conds, xvals, zvals, new_times = get_vals_conds(   q_reaper, m, qfs, qs, qfname, fcn,
                                                                axes, x_axis, dsample, z_axis)

    if x_axis == 'time':
        xvals = np.linspace(0, dsample-1, dsample)

    for i in np.ndindex(fd_vals.shape[:-2]):
        savename = 'model_{}_{}_{}_{}_vs-{}-and-{}-3D_{}.png'.format(m, process, t_or_w, savelab, x_axis, z_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) for n, cond in enumerate(conds)]))
        print(savename)

        X, Y = np.meshgrid(xvals,zvals)

        elev, azim = 20, 60
        ax = plt.axes(projection='3d')
        ax.plot_surface(X, Y, fd_vals[i], rstride=1, cstride=1, cmap=cmap, edgecolor='none')
        ax.set_xlabel(xlab)
        ax.set_ylabel(zlab)
        ax.set_zlabel(qlab)
        ax.view_init(elev, azim)

        plt.tight_layout()
        plt.savefig(savename)
        plt.close()

# plot_fcn_of_quants_3D(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_tot'], cov,
#                         base_axes, x_axis='To', dsample=None, z_axis='L',
#                         xlab=r'$T^o_2$',
#                         qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                         zlab='L', cmap='viridis', t_or_w='w',
#                         logscale=0, savelab='ifio-mitot-corr')

# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'], cov,
#                         time_axes, x_axis='time', dsample=51, z_axis='To',
#                         xlab=r'$50*t/(t^{NESS})$',
#                         qlab=r'$\chi(I^{i \to o}, I^{tot}(i;o))$',
#                         zlab='To', cmap='viridis', t_or_w='t',
#                         logscale=0, savelab='ifio-mitot-corr')






'''
################################################################################
S6, Static plot histogram(quant) vs param, 3D
################################################################################
and can also be used for
################################################################################
A3, Animation of histogram(quant) vs param/time animated in param/time, 3D
################################################################################
'''
qcmaps = {  'if_io':    plt.cm.cool,
            'if_oi':    plt.cm.PuRd,
            'mi_inst':  plt.cm.YlGnBu,
            'mi_tot':   plt.cm.PuBuGn,
            'ep':       plt.cm.YlOrRd,
            'q':        plt.cm.Wistia,
            's':        plt.cm.summer}

def get_cmap (base, data):
    return base(plt.Normalize(0,np.max(data))(data))

def plot_3d_dists(  q_reaper, m, qfs, q, qfname,
                    axes, x_axis, xlab, dsample,
                    nbins, c, qlab, t_or_w,
                    anim_axis=None, zlim=200, process=''):
    '''
    q_reaper (fcn): to use to extract quants,
    m (str): model name,
    qfs (nparray): storing names of files with desired quantities - shape dims ,
    q (str): name of quantity to extract,
    axes (dict): of 'this': axis along which 'this' varies,
    x_axis (str): in 3D histogram,
    dsample (float): if provided, how much to downsample anim_axis,
    nbins (int): number of bins for histogram,
    c (str): color for histogram bars,
    qlab (str): label for x-axis of histogram,
    t_or_w (str): 't' if for trajectories, 'w' if for windows,
    anim_axis (str): along which axis to animate,
    '''
    if x_axis == 'time':
        vals = get_quants(q_reaper, m, qfs, q, qfname, dsample, axes, x_axis)
    else:
        vals = get_quants(q_reaper, m, qfs, q, qfname, None, None, None)
    vals, new_axes = roll_vals_axes(vals, axes, x_axis, axes['num_axes']-1)

    if anim_axis is None:
        qmin, qmax = np.amin(vals), np.amax(vals)
        bin_list = np.linspace(qmin, qmax, nbins + 1)
    else:
        vals, new_axes = roll_vals_axes(vals, new_axes, anim_axis, axes['num_axes']-2)


    new_axes = dict(sorted(new_axes.items(), key = lambda i: i[1]))
    conds = np.array([i for i in new_axes.keys()])
    conds = np.delete(conds, [new_axes[x_axis], new_axes['traj'], new_axes['num_axes']])

    if 'time' in new_axes and x_axis != 'time':
        vals = downsample_along_axis(vals, dsample, new_axes['time'])
        tvals = np.linspace(0, dsample-1, dsample)

    if x_axis != 'time':
        xvals = conds_dict[m][x_axis]
    else:
        xvals = np.linspace(0, dsample-1, dsample)

    for i in np.ndindex(vals.shape[:-2]):
        savename = 'model_{}_{}_{}_3Dhist-{}_vs-{}_for-{}.png'.format(m, process, t_or_w, q, x_axis, '_'.join(['{}-{}'.format(cond, conds_dict[m][cond][i[n]]) if cond in conds_dict[m] else '{}-{}'.format('time', tvals[i[n]]) for n, cond in enumerate(conds)]))
        print(savename)

        if anim_axis is not None:
            qmin, qmax = np.amin(vals[i[:-1]]), np.amax(vals[i[:-1]])
            bin_list = np.linspace(qmin, qmax, nbins + 1)

        hists = np.zeros((vals[i].shape[0], nbins))
        for j in range(len(xvals)):
            hist, _ = np.histogram(vals[i][j], bin_list)
            hists[j] = hist

        x = 0.5*(bin_list[:-1] + bin_list[1:]) # centers of bins
        y = xvals
        X, Y = np.meshgrid(x,y)
        ax = plt.axes(projection='3d')
        x_data = X.flatten()
        y_data = Y.flatten()
        z_data = hists.flatten()

        color_data = z_data
        cmap = get_cmap(qcmaps[q], color_data)
        # cmap = plt.cm.cool(plt.Normalize(0,np.max(z_data))(z_data))
        ax.bar3d(   x_data, y_data, np.zeros(len(z_data)), x[1]-x[0],
                    y[1]-y[0], z_data, color = cmap, alpha=0.9)
        ax.set_zlim(top=zlim)
        ax.set_xlabel(qlab)
        ax.set_ylabel(xlab)
        plt.tight_layout()
        elev, azim = 20, 60
        ax.view_init(elev, azim)
        plt.savefig(savename)
        plt.close()

# plot_3d_dists(      get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                     timextraj_axes, 'time', xlab=r'$t$', dsample=51,
#                     nbins=50, c=qcolors['if_io'],
#                     qlab=r'$I^{i \to o}$', t_or_w='t',
#                     anim_axis='To', zlim=10)

# plot_3d_dists(      get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                     timextraj_axes, 'To', xlab=r'$T^o_2$', dsample=51,
#                     nbins=50, c=qcolors['if_io'],
#                     qlab=r'$I^{i \to o}$', t_or_w='t',
#                     anim_axis='time', zlim=10)

# plot_3d_dists(      get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                     traj_axes, 'To', xlab=r'$T^o_2$',
#                     dsample=None, nbins=50, c=qcolors['ep'],
#                     qlab=r'$<\sigma>$', t_or_w='t')
#
# plot_3d_dists(      get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                     traj_axes, 'To', xlab=r'$T^o_2$',
#                     dsample=None, nbins=50, c=qcolors['ep'],
#                     qlab=r'$<\sigma>$', t_or_w='w')


# use glob to generate gifs,
# need to rewrite these functions

'''Converting .pngs from S6, Relaxation --> .gif of A3, Relaxation'''
def get_img_filenames(modelname, quant, param, id, num_msg, anim_dim, num):
    files = []
    for i in range(num):
        if anim_dim == 'param':
            files.append('{}_{}-{}-vs-{}{}-time-3D_for_{}-{}_{}-{}.png'.format(modelname, quant, 'hist', param, str(i).zfill(4), 'id', id, 'L', num_msg))
        elif anim_dim == 'time':
            files.append('{}_{}-{}-vs-{}-time{}-3D_for_{}-{}_{}-{}.png'.format(modelname, quant, 'hist', param, str(i).zfill(4), 'id', id, 'L', num_msg))
    return files


def convert_all_togifs(modelname, quants, param, ids, num_msgs, anim_dim, num, fps=5):
    for quant, id, num_msg in product(quants, ids, num_msgs):
        images = []
        filenames = get_img_filenames(modelname, quant, param, id, num_msg, anim_dim, num)
        for filename in filenames:
            images.append(imageio.imread(filename))
        if anim_dim == 'time':
            imageio.mimsave('{}_{}-{}-vs-{}-timeanim-3D_for_{}-{}_{}-{}.gif'.format(modelname, quant, 'hist', param, 'id', id, 'L', num_msg), images, fps=fps)
        if anim_dim == 'param':
            imageio.mimsave('{}_{}-{}-vs-{}anim-time-3D_for_{}-{}_{}-{}.gif'.format(modelname, quant, 'hist', param, 'id', id, 'L', num_msg), images, fps=fps)


# convert_all_togifs('model_r1', ['ep', 'if_io', 'mi'], 'eps', R1_initdists, n_msgs[:2], 'param', 11)
# convert_all_togifs('model_r1', ['ep', 'if_io', 'mi'], 'To', R1_initdists, n_msgs[:2], 'param', 11)

# convert_all_togifs('model_r1', ['ep', 'if_io', 'mi'], 'eps', R1_initdists, n_msgs, 'time', 51)
# convert_all_togifs('model_r1', ['ep', 'if_io', 'mi'], 'To', R1_initdists, n_msgs, 'time', 51)
