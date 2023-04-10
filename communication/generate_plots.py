from plots import *


'''
################################################################################
S1, Cness vs params, 2D
################################################################################
'''
plt.rc('legend',fontsize='x-small')
# '''
# --------------------------------------------------------------------------------
# Cness vs To, with markers as H and colors as L
# --------------------------------------------------------------------------------
# '''
# plot_Cness_vparams( 'nr', x_axis='To', m_axis='H', c_axis='L',
#                     ms=markers, cs=pcolors['L'],
#                     xname=r'$T^o_2$', mname=r'$H$', cname=r'$L$',
#                     t_or_w='t', Cnesslog=1, legloc = 'upper right')
#
# plot_Cness_vparams( 'r1', x_axis='To', m_axis='H', c_axis='L',
#                     ms=markers, cs=pcolors['L'],
#                     xname=r'$T^o_2$', mname=r'$H$', cname=r'$L$',
#                     t_or_w='t', Cnesslog=1, legloc=None)
# '''
# --------------------------------------------------------------------------------
# Cness vs sr, with markers as H and colors as L
# --------------------------------------------------------------------------------
# '''
# plot_Cness_vparams( 'nr', x_axis='sr', m_axis='H', c_axis='L',
#                     ms=markers, cs=pcolors['L'],
#                     xname=r'$f_s$', mname=r'$H$', cname=r'$L$',
#                     t_or_w='t', Cnesslog=1)


# '''
# --------------------------------------------------------------------------------
# Cness vs L, with markers as H and colors as sr
# --------------------------------------------------------------------------------
# '''
# plot_Cness_vparams( 'nr', x_axis='L', m_axis='H', c_axis='sr',
#                     ms=markers, cs=pcolors['sr_nr'],
#                     xname=r'$L$', mname=r'$H$', cname=r'$f_s$',
#                     t_or_w='t', Cnesslog=0)
#
# plot_Cness_vparams( 'r1', x_axis='L', m_axis='H', c_axis='sr',
#                     ms=markers, cs=pcolors['sr_r1'],
#                     xname=r'$L$', mname=r'$H$', cname=r'$f_s$',
#                     t_or_w='t', Cnesslog=0)
# '''
# --------------------------------------------------------------------------------
# Cness vs L, with markers as H and colors as To
# --------------------------------------------------------------------------------
# '''
# plot_Cness_vparams( 'nr', x_axis='L', m_axis='H', c_axis='To',
#                     ms=markers, cs=pcolors['To'],
#                     xname=r'$L$', mname=r'$H$', cname=r'$T^o_2$',
#                     t_or_w='t', Cnesslog=0)
#
# plot_Cness_vparams( 'r1', x_axis='L', m_axis='H', c_axis='To',
#                     ms=markers, cs=pcolors['To'],
#                     xname=r'$L$', mname=r'$H$', cname=r'$T^o_2$',
#                     t_or_w='t', Cnesslog=1, legloc=None)










'''
################################################################################
A1, Animated histograms, 2D
################################################################################
'''


# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each window, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma_w>$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each window, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}_w$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each window, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}_w(i;o)$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each window, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}_w(i;o)$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each window, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma_w>$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each window, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}_w$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each window, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}_w(i;o)$', t_or_w='w', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each window, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot', 'nr_wind_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}_w(i;o)$', t_or_w='w', process = 'ness')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each trajectory, NESS, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep', 'r1_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(10)>$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each trajectory, NESS, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io', 'r1_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each trajectory, NESS, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst', 'r1_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each trajectory, NESS, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot', 'r1_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;10)$', t_or_w='t', process = 'ness')
#
#
#
#
#
#
# #
# #
# #
# #
# #
# #
# #
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each trajectory, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(10)>$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each trajectory, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each trajectory, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each trajectory, NESS, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each trajectory, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(10)>$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each trajectory, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each trajectory, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;10)$', t_or_w='t', process = 'ness')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each trajectory, NESS, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot', 'nr_traj_ness_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;10)$', t_or_w='t', process = 'ness')
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
#
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each relaxation, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep', 'r1_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(t_r)>$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each relaxation, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io', 'r1_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each relaxation, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst', 'r1_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each relaxation, r1, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot', 'r1_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each relaxation, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(t_r)>$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each relaxation, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each relaxation, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each relaxation, nr, anim with To
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='To',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP by the end of each relaxation, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['ep'], fps=5,
#                     qlab=r'$<\sigma(t_r)>$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io by the end of each relaxation, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['if_io'], fps=5,
#                     qlab=r'$I^{i \to o}(t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst by the end of each relaxation, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_inst'], fps=5,
#                     qlab=r'$I^{inst}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot by the end of each relaxation, nr, anim with sr
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot', 'nr_fin_rlx_quants',
#                     traj_axes, anim_axis='sr',
#                     dsample=None, nbins=50, c=qcolors['mi_tot'], fps=5,
#                     qlab=r'$I^{tot}(i;o;t_r)$', t_or_w='t', process = 'fin_rlx')











# '''
# --------------------------------------------------------------------------------
# Distribution of total EP while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep', 'r1_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['ep'], fps=10,
#                     qlab=r'$<\sigma(t_r)>$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io', 'r1_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['if_io'], fps=10,
#                     qlab=r'$I^{i \to o}(t_r)$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst', 'r1_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['mi_inst'], fps=10,
#                     qlab=r'$I^{inst}(i;o;t_r)$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot', 'r1_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['mi_tot'], fps=10,
#                     qlab=r'$I^{tot}(i;o;t_r)$', t_or_w='t', process='rlx')
#
#
#
#
#
#
#
#
#
# '''
# --------------------------------------------------------------------------------
# Distribution of total EP while relaxing, nr, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep', 'nr_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['ep'], fps=10,
#                     qlab=r'$<\sigma(t_r)>$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total IF_io while relaxing, nr, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io', 'nr_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['if_io'], fps=10,
#                     qlab=r'$I^{i \to o}(t_r)$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_inst while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst', 'nr_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['mi_inst'], fps=10,
#                     qlab=r'$I^{inst}(i;o;t_r)$', t_or_w='t', process='rlx')
# '''
# --------------------------------------------------------------------------------
# Distribution of total MI_tot while relaxing, r1, anim with time
# --------------------------------------------------------------------------------
# '''
# animate_2d_dists(get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot', 'nr_rlx_quants',
#                     timextraj_axes, anim_axis='time',
#                     dsample=51, nbins=50, c=qcolors['mi_tot'], fps=10,
#                     qlab=r'$I^{tot}(i;o;t_r)$', t_or_w='t', process='rlx')










'''
################################################################################
S2, Static plot of f(quants) vs param, 2D
################################################################################
'''



# '''
# --------------------------------------------------------------------------------
# Precision of quantity by end of NESS windows, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(\sigma_w)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(\sigma_w)$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{i \to o}_w)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{i \to o}_w)$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{inst}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{inst}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{tot}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{tot}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='precision-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Mean of quantity by end of NESS windows, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<\sigma_w>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<\sigma_w>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{i \to o}_w>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{i \to o}_w>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{inst}_w(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{inst}_w(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{tot}_w(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{tot}_w(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='avg-mi_tot', process='ness')
#
#
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS windows, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
#

plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', nr_wind_ness_qfs,
                                ['ep', 'mi_tot'], 'r1_traj_ness_quants', squared_pearson_cc,
                                base_axes, m_axis='L', c_axis='To',
                                ms=markers, cs=pcolors['To'],
                                qlab=r'$\chi^2(\sigma_w, I^{tot}_w(i;o))$',
                                mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
                                ylog=0, savelab='pcc-ep-mi_tot', process='ness')

plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
                                ['ep', 'mi_tot'], 'nr_wind_ness_quants', squared_pearson_cc,
                                base_axes, m_axis='L', c_axis='sr',
                                ms=markers, cs=pcolors['sr_nr'],
                                qlab=r'$\chi^2(\sigma_w, I^{tot}_w(i;o))$',
                                mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
                                ylog=0, savelab='pcc-ep-mi_tot', process='ness')

plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
                                ['ep', 'mi_tot'], 'nr_wind_ness_quants', squared_pearson_cc,
                                base_axes, m_axis='L', c_axis='To',
                                ms=markers, cs=pcolors['To'],
                                qlab=r'$\chi^2(\sigma_w, I^{tot}_w(i;o))$',
                                mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
                                ylog=0, savelab='pcc-ep-mi_tot', process='ness')

# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_wind_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}_w, I^{inst}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_wind_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}_w, I^{inst}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_wind_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}_w, I^{tot}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_wind_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}_w, I^{tot}_w(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# IFT test of quantity by end of NESS windows, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-\sigma_w}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'ep', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-\sigma_w}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{i \to o}_w}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'if_io', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{i \to o}_w}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-if_io', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{inst}_w(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_inst', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{inst}_w(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{tot}_w(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                 'mi_tot', 'nr_wind_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{tot}_w(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 'w',
#                                 ylog=0, savelab='ift-mi_tot', process='ness')


# '''
# --------------------------------------------------------------------------------
# Precision of quantity by end of NESS trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='ness')
#

# '''
# --------------------------------------------------------------------------------
# Mean of quantity by end of NESS trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='ness')
#
#
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# IFT test of quantity by end of NESS trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'ep', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='ness')

# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'if_io', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_inst', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='ness')
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                 'mi_tot', 'nr_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='ness')


# '''
# --------------------------------------------------------------------------------
# Precision of quantity by end of rlx trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='fin_rlx')
#
#
# '''
# --------------------------------------------------------------------------------
# Mean of quantity by end of rlx trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='fin_rlx')
#
#
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 ['if_io', 'mi_inst'], 'nr_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 ['if_io', 'mi_tot'], 'nr_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# IFT test of quantity by end of rlx trajectories, nr, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'ep', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'if_io', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_inst', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                 'mi_tot', 'nr_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='sr',
#                                 ms=markers, cs=pcolors['sr_nr'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$f_s$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='fin_rlx')


# '''
# --------------------------------------------------------------------------------
# Precision of quantity by end of NESS trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'ep', 'r1_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'if_io', 'r1_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_inst', 'r1_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_tot', 'r1_traj_ness_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Mean of quantity by end of NESS trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'ep', 'r1_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'if_io', 'r1_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_inst', 'r1_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_tot', 'r1_traj_ness_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='ness')
#
#
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 ['if_io', 'mi_inst'], 'r1_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 ['if_io', 'mi_tot'], 'r1_traj_ness_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# IFT test of quantity by end of NESS trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'ep', 'r1_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'if_io', 'r1_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_inst', 'r1_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='ness')
#
#
# plot_2d_fcn_of_quants_vCness(   get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                 'mi_tot', 'r1_traj_ness_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Precision of quantity by end of rlx trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'ep', 'r1_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(\sigma)$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'if_io', 'r1_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{i \to o})$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_inst', 'r1_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_tot', 'r1_fin_rlx_quants', precision,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$prec(I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='precision-mi_tot', process='fin_rlx')
#
#
# '''
# --------------------------------------------------------------------------------
# Mean of quantity by end of rlx trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'ep', 'r1_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<\sigma>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'if_io', 'r1_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{i \to o}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_inst', 'r1_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{inst}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_tot', 'r1_fin_rlx_quants', avg,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<I^{tot}(i;o)>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='avg-mi_tot', process='fin_rlx')
#
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 ['if_io', 'mi_inst'], 'r1_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 ['if_io', 'mi_tot'], 'r1_fin_rlx_quants', squared_pearson_cc,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='pcc-if_io-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# IFT test of quantity by end of rlx trajectories, r1, marker as L, color as T
# --------------------------------------------------------------------------------
# '''
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'ep', 'r1_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-\sigma}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-ep', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'if_io', 'r1_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{i \to o}}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-if_io', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_inst', 'r1_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{inst}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_inst', process='fin_rlx')
#
# plot_2d_fcn_of_quants_vCness(   get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs,
#                                 'mi_tot', 'r1_fin_rlx_quants', ift,
#                                 base_axes, m_axis='L', c_axis='To',
#                                 ms=markers, cs=pcolors['To'],
#                                 qlab=r'$<e^{-I^{tot}(i;o)}>$',
#                                 mname=r'$L$', cname=r'$T^o_2$', t_or_w = 't',
#                                 ylog=0, savelab='ift-mi_tot', process='fin_rlx')


'''
################################################################################
S3, Static plot of multiple f(quants) vs param, 2D
################################################################################
'''
# '''
# --------------------------------------------------------------------------------
# TUR test of quantity by end of NESS windows, nr, windows
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['if_io', 'ep'],
#                                     'nr_wind_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$1/prec(I^{i \to o}_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$f_s$', t_or_w = 'w',
#                                     ylog=0, savelab='tur-if_io', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['if_io', 'ep'],
#                                     'nr_wind_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(I^{i \to o}_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$T^o_2$', t_or_w = 'w',
#                                     ylog=0, savelab='tur-if_io', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'ep'],
#                                     'nr_wind_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$1/prec(\sigma_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$f_s$', t_or_w = 'w',
#                                     ylog=0, savelab='tur-ep', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'ep'],
#                                     'nr_wind_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(\sigma_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$T^o_2$', t_or_w = 'w',
#                                     ylog=0, savelab='tur-ep', process='ness')
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS windows, nr, windows
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_wind_ness_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$f_s$', t_or_w = 'w',
#                                     ylog=0, savelab='avg-quants', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_wind_ness_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$T^o_2$', t_or_w = 'w',
#                                     ylog=0, savelab='avg-quants', process='ness')
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS windows, nr, windows
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_wind_ness_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$f_s$', t_or_w = 'w',
#                                     ylog=0, savelab='precision-quants', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_wind_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_wind_ness_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$T^o_2$', t_or_w = 'w',
#                                     ylog=0, savelab='precision-quants', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# TUR test of quantity by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['if_io', 'ep'],
#                                     'nr_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$1/prec(I^{i \to o}_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='tur-if_io', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['if_io', 'ep'],
#                                     'nr_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(I^{i \to o}_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='tur-if_io', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'ep'],
#                                     'nr_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$1/prec(\sigma_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='tur-ep', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'ep'],
#                                     'nr_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(\sigma_w)$',r'$2/<\sigma_w>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='tur-ep', process='ness')
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_traj_ness_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_traj_ness_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='ness')
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_traj_ness_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='ness')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'nr', nr_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_traj_ness_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='ness')



# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_fin_rlx_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='fin_rlx')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_fin_rlx_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$<\sigma_w>$', r'$<I^{i \to o}_w>$', r'$<I^{inst}_w(i;o)>$', r'$<I^{tot}_w(i;o)>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_fin_rlx_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='sr',
#                                     ms=qmarkers, cs=pcolors['sr_nr'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$f_s$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='fin_rlx')
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'nr_fin_rlx_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$prec(\sigma_w)$', r'$prec(I^{i \to o}_w)$', r'$prec(I^{inst}_w(i;o))$', r'$prec(I^{tot}_w(i;o))$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='fin_rlx')


# '''
# --------------------------------------------------------------------------------
# TUR test of quantity by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                     ['if_io', 'ep'],
#                                     'r1_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(I^{i \to o})$',r'$2/<\sigma>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='tur-if_io', process='ness')
#
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                     ['ep', 'ep'],
#                                     'r1_traj_ness_quants',
#                                     [lack_of_precision, twice_reciprocal_avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$1/prec(\sigma)$',r'$2/<\sigma>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='tur-ep', process='ness')
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'r1_traj_ness_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$<\sigma>$', r'$<I^{i \to o}>$', r'$<I^{inst}(i;o)>$', r'$<I^{tot}(i;o)>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='ness')
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_ness_quants, 'r1', r1_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'r1_traj_ness_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$prec(\sigma)$', r'$prec(I^{i \to o})$', r'$prec(I^{inst}(i;o))$', r'$prec(I^{tot}(i;o))$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='ness')
#
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'r1', r1_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'r1_fin_rlx_quants',
#                                     [avg, avg, avg, avg],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$<\sigma>$', r'$<I^{i \to o}>$', r'$<I^{inst}(i;o)>$', r'$<I^{tot}(i;o)>$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='avg-quants', process='fin_rlx')
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
#
# plot_2d_mult_fcn_of_quants_vCness(  get_all_fin_rlx_quants, 'r1', r1_traj_ness_qfs,
#                                     ['ep', 'if_io', 'mi_inst', 'mi_tot'],
#                                     'r1_fin_rlx_quants',
#                                     [precision, precision, precision, precision],
#                                     base_axes, c_axis='To',
#                                     ms=qmarkers, cs=pcolors['To'],
#                                     fqlabs=[r'$prec(\sigma)$', r'$prec(I^{i \to o})$', r'$prec(I^{inst}(i;o))$', r'$prec(I^{tot}(i;o))$'],
#                                     cname=r'$T^o_2$', t_or_w = 't',
#                                     ylog=0, savelab='precision-quants', process='fin_rlx')


'''
################################################################################
S4, Static plot heatmapping f(quants) vs param and time, 2D
################################################################################
'''

# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'if_io',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_inst',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_tot',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'if_io',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_inst',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_tot',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_inst'],
#                             'nr_wind_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='w',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_tot'],
#                             'nr_wind_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='w',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='ness')


# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, ['if_io', 'mi_inst'],
#                             'nr_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'nr', nr_traj_ness_qfs, ['if_io', 'mi_tot'],
#                             'nr_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='ness')


# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='fin_rlx')


# '''
# --------------------------------------------------------------------------------
# Averages of quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')

# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')

# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')

# '''
# --------------------------------------------------------------------------------
# Precisions of quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')

# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, ['if_io', 'mi_inst'],
#                             'r1_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_heatmap( get_all_ness_quants, 'r1', r1_traj_ness_qfs, ['if_io', 'mi_tot'],
#                             'r1_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'r1_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_heatmap( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'r1_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='fin_rlx')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'r1_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_heatmap( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'r1_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')

'''
################################################################################
S5, Static plot surface of f(quants) vs param and time, 3D
################################################################################
'''

# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'if_io',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_inst',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_tot',
#                             'nr_wind_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='w',
#                             logscale=0, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'if_io',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_inst',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_tot',
#                             'nr_wind_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='w',
#                             logscale=0, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS windows, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_inst'],
#                             'nr_wind_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='w',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_wind_ness_qfs, ['if_io', 'mi_tot'],
#                             'nr_wind_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='w',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot',
#                             'nr_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot',
#                             'nr_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, ['if_io', 'mi_inst'],
#                             'nr_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'nr', nr_traj_ness_qfs, ['if_io', 'mi_tot'],
#                             'nr_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='sr',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='fin_rlx')


# '''
# --------------------------------------------------------------------------------
# Averages of quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$f_s$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                             'nr_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')

# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities during rlx trajectories, nr
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='sr',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$f_s$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'nr_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot',
#                             'r1_traj_ness_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot',
#                             'r1_traj_ness_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='ness')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of NESS trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, ['if_io', 'mi_inst'],
#                             'r1_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='ness')
#
# plot_fcn_of_quants_3D( get_all_ness_quants, 'r1', r1_traj_ness_qfs, ['if_io', 'mi_tot'],
#                             'r1_traj_ness_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='ness')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-ep', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_fin_rlx_quants', avg,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$L$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=1, savelab='avg-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-ep', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-if_io', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_fin_rlx_quants', precision,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=1, savelab='precision-mi_tot', process='fin_rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities by end of rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'r1_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_inst', process='fin_rlx')
#
# plot_fcn_of_quants_3D( get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'r1_fin_rlx_quants', squared_pearson_cc,
#                             base_axes, x_axis='To', dsample=None, z_axis='L',
#                             xlab=r'$T^o_2$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$L$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=1, savelab='pcc-if_io-mi_tot', process='fin_rlx')
#
#
# '''
# --------------------------------------------------------------------------------
# Averages of quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<\sigma_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{i \to o}_w>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{inst}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_rlx_quants', avg,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$<I^{tot}_w(i;o)>$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['mean'], t_or_w='t',
#                             logscale=0, savelab='avg-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Precisions of quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(\sigma_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-ep', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{i \to o}_w)$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-if_io', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{inst}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                             'r1_rlx_quants', precision,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$prec(I^{tot}_w(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['precision'], t_or_w='t',
#                             logscale=0, savelab='precision-mi_tot', process='rlx')
#
# '''
# --------------------------------------------------------------------------------
# Pearson correlation coefficient between two quantities during rlx trajectories, r1
# --------------------------------------------------------------------------------
# '''
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_inst'],
#                             'r1_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{inst}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_inst', process='rlx')
#
# plot_fcn_of_quants_3D( get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, ['if_io', 'mi_tot'],
#                             'r1_rlx_quants', squared_pearson_cc,
#                             time_axes, x_axis='time', dsample=51, z_axis='To',
#                             xlab=r'$t$',
#                             qlab=r'$\chi^2(I^{i \to o}, I^{tot}(i;o))$',
#                             zlab=r'$T^o_2$', cmap=fcmaps['pearson_cc'], t_or_w='t',
#                             logscale=0, savelab='pcc-if_io-mi_tot', process='rlx')

'''
################################################################################
S6, Static plot histogram(quant) vs param, 3D
################################################################################
and can also be used for
################################################################################
A3, Animation of histogram(quant) vs param/time animated in param/time, 3D
################################################################################
'''

# plot_3d_dists(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'ep',
#                 'nr_wind_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='w',
#                 anim_axis='sr', zlim=200, process='ness')

# plot_3d_dists(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'if_io',
#                 'nr_wind_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='w',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_inst',
#                 'nr_wind_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='w',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'nr', nr_wind_ness_qfs, 'mi_tot',
#                 'nr_wind_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='w',
#                 anim_axis='sr', zlim=200, process='ness')


# plot_3d_dists(  get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'ep',
#                 'nr_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'if_io',
#                 'nr_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_inst',
#                 'nr_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'nr', nr_traj_ness_qfs, 'mi_tot',
#                 'nr_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                 'nr_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                 'nr_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                 'nr_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                 'nr_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'ep',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'if_io',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_inst',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'nr', nr_traj_rlx_qfs, 'mi_tot',
#                 'nr_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='sr', zlim=200, process='rlx')


# plot_3d_dists(  get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'ep',
#                 'r1_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'if_io',
#                 'r1_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_inst',
#                 'r1_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_ness_quants, 'r1', r1_traj_ness_qfs, 'mi_tot',
#                 'r1_traj_ness_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='ness')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                 'r1_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'r1', nr_traj_rlx_qfs, 'if_io',
#                 'r1_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                 'r1_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='fin_rlx')
#
# plot_3d_dists(  get_all_fin_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                 'r1_fin_rlx_quants', traj_axes, 'To', xlab=r'$T^o_2$',
#                 dsample=None, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='L', zlim=200, process='fin_rlx')

# plot_3d_dists(  get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'ep',
#                 'r1_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['ep'],
#                 qlab=r'$<\sigma>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'if_io',
#                 'r1_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['if_io'],
#                 qlab=r'$<I^{i \to o}>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_inst',
#                 'r1_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_inst'],
#                 qlab=r'$<I^{inst}(i;o)>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
#
# plot_3d_dists(  get_all_rlx_quants, 'r1', r1_traj_rlx_qfs, 'mi_tot',
#                 'r1_rlx_quants', timextraj_axes, 'time', xlab=r'$t$',
#                 dsample=51, nbins=50, c=qcolors['mi_tot'],
#                 qlab=r'$<I^{tot}(i;o)>$', t_or_w='t',
#                 anim_axis='To', zlim=200, process='rlx')
