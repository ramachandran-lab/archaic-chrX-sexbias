#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 3.

    The effects of dominance and sex-bias on chromosome X archaic coverage.
    A, B
        Mean coverage ratios between a simulated autosome and chromosome X,  estimated from 10,000 pairwise combinations of simulated autosomal and X chromosomes.
        Filled points reflect additive variants ($h=0.5$); hollow points reflect recessive variants ($h=0.0$).
        Horizontal lines denote autosome to chromosome X ratio of archaic haplotypes in initial introgression pulse. See Eqn ref{eqn:haplotypes}.
    C, D
        Filled violins depict additive variants; hollow violins depict recessive variants.  Horizontal lines denote magnitude of initial introgression pulse.
"""

import extract_pI_data as pI
import extract_pI_data_sexbias as sexpI

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as colors
import numpy as np
import pandas as pd


# %% Import data

path_to_pI_files = 'fig3_data/nonsex/'
pI_df = pI.compile_df_from_pI_files(path_to_pI_files)

path_to_sexpI_files = 'fig3_data/sex/'
sexpI_df = sexpI.compile_df_from_pI_files(path_to_sexpI_files)

megadf = pd.concat([pI_df, sexpI_df])
megadf['target_male_frac'] = megadf['target_male_frac'].fillna(0.5)


chr1_add_df = sexpI.subset_pI(megadf, c='1', d='0.5SLiM', s=1, n=0.01)
chr1_rec_df = sexpI.subset_pI(megadf, c='1', d='0.0', s=1, n=0.01)
chrX_add_df = sexpI.subset_pI(megadf, c='X', d='0.5SLiM', s=1, n=0.01)
chrX_rec_df = sexpI.subset_pI(megadf, c='X', d='0.0', s=1, n=0.01)


# %% Calculate aut:chrX ratio from sim data

def bs_autx_prop(aut_pIs, x_pIs, ratio_handling='ratioOfMeans', n=1000, B=2000):
    """
    Bootstrap the aut:chrX coverage proportion from simulation data.

    Parameters
    ----------
    aut_pIs : np array
        List of simulated coverages on autosomal haplotypes.
    x_pIs : np array
        List of simulated coverages on chromosome X haplotypes.
    n : int, optional
        Number of ratios to generate. The default is 1000.
    B : int, optional
        Number of bootstrap replications. The default is 2000.

    Returns
    -------
    this_mean_stat : float
        Mean of all bootstrapped aut:chrX ratios.
    float
        Standard error of bootstrapped aut:chrX ratios.
    tuple
        (B, (len(aut_pIs), len(x_pIs))).
    these_cis : tuple
        2.5 and 97.5%th percentiles of bootstrapped aut:chrX ratios.

    """
    if n is None:
        n = min(len(aut_pIs), len(x_pIs))
        print(n)
    assert min(len(aut_pIs), len(x_pIs)) >= n

    a_choices = np.random.choice(aut_pIs, size=(B, n), replace=True)
    x_choices = np.random.choice(x_pIs, size=(B, n), replace=True)

    if ratio_handling == 'ratioOfMeans':
        a_means = np.mean(a_choices, axis=1)
        x_means = np.mean(x_choices, axis=1)
        bs_ratios = a_means / x_means
    elif ratio_handling == 'meanOfRatios':
        ratios_for_bs = a_choices / x_choices
        bs_ratios = np.mean(ratios_for_bs, axis=1)
    this_mean_stat = np.mean(bs_ratios)
    these_cis = np.percentile(bs_ratios, (2.5, 97.5))

    return this_mean_stat, np.std(bs_ratios), (len(bs_ratios), (len(aut_pIs), len(x_pIs))), these_cis

mfrac_list = (1, 0.5, 0)

n = None  # number of aut:chrX ratios to generate
additive_ratio_stats = {}
recessive_ratio_stats = {}
for m in mfrac_list:
    aadf = sexpI.subset_pI(chr1_add_df, m=m)
    xadf = sexpI.subset_pI(chrX_add_df, m=m)
    ardf = sexpI.subset_pI(chr1_rec_df, m=m)
    xrdf = sexpI.subset_pI(chrX_rec_df, m=m)

    astats = bs_autx_prop(np.asarray(aadf['pI']), np.asarray(xadf['pI']), n=n)
    additive_ratio_stats.update({m: astats})

    rstats = bs_autx_prop(np.asarray(ardf['pI']), np.asarray(xrdf['pI']), n=n)
    recessive_ratio_stats.update({m: rstats})
    print(f"Finished mfrac {m}.")


# %% Plotting preparation

cmap = colors.get_cmap(name='tab20')
blue = cmap(0)
orange = cmap(2)


def pI_by_targetmfrac_violin(ax, input_df,
                             mfrac_list=[0, 0.2, 0.4, 0.6, 0.8, 1],
                             width=0.2, face_edge=None, return_data=False,
                             flip=False):
    list_of_box_data = []
    for m in mfrac_list:
        df = sexpI.subset_pI(input_df, m=m)
        these_pIs = df['pI'].tolist()
        list_of_box_data.append(these_pIs)

    # Recolor violins for dominance and chromosome type
    if face_edge is not None:
        face, edge = face_edge
        parts = ax.violinplot(list_of_box_data, positions=mfrac_list,
                              widths=width, showmeans=False, points=300)
        for pc in parts['bodies']:
            pc.set_facecolor(face)
            pc.set_edgecolor(edge)
        # https://stackoverflow.com/questions/26291479/changing-the-color-of-matplotlibs-violin-plots#comment95989703_26291582
        parts['cbars'].set_edgecolor(edge)
        parts['cmaxes'].set_color(edge)
        parts['cmins'].set_edgecolors(edge)
    else:
        ax.violinplot(list_of_box_data, positions=mfrac_list, widths=width,
                      showmeans=False, points=300)
    ax.set_xticks(mfrac_list)
    if flip:
        ax.set_xticklabels(['male\nintrogressors\n$p = 0$',
                            'no\nsexbias\n$p = 0.5$',
                            'female\nintrogressors\n$p = 1$'],
                           fontsize=8)
    # # initial introg line(s)
    # ax.axhline(y=0.05, linestyle='--', color='r')
    if return_data:
        return ax, list_of_box_data
    else:
        return ax


# %% * Plot Fig 3

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10


def plot_sexbias_violins_with_ratio(add_stats, rec_stats, mfracs=(1, 0.5, 0),
                                    errs='e'):
    """
    Generate Fig 3.

    Four-panel figure with panel of coverage violins and panel of
    mean aut:chrX ratios with errorbars for additive simulations,
    and same for recessive simulations.


    Parameters
    ----------
    add_stats : dict
        Coverage statistics from additive simulations, calculated as above.
    rec_stats : dict
        Coverage statistics from recessive simulations, calculated as above.
    mfracs : tuple of floats, optional
        Male fractions, left to right. The default is (1, 0.5, 0).
    errs : char, optional
        Type of errorbar for aut:chrX ratio.  'e' for empirical 95% CI from
        bootstraps stored in add/rec_stats dictionary.  'n' for 95% CI from
        normal distribution with standard error stored in add/rec_stats dictionary.
        The default is 'e'.

    Returns
    -------
    N/A
    """
    fig, axes = plt.subplots(2, 2, gridspec_kw={'height_ratios': [1, 3]},
                             sharex='col', sharey='row', figsize=(5.2, 5.2))
    # https://stackoverflow.com/a/63770891/18621926
    raxadd = axes[0][0]
    raxrec = axes[0][1]
    axadd = axes[1][0]
    axrec = axes[1][1]

    # violins
    axadd = pI_by_targetmfrac_violin(axadd, chr1_add_df,
                                     face_edge=(blue, blue),
                                     mfrac_list=mfracs, flip=True)
    axadd = pI_by_targetmfrac_violin(axadd, chrX_add_df,
                                     face_edge=(orange, orange),
                                     mfrac_list=mfracs, flip=True)
    axrec = pI_by_targetmfrac_violin(axrec, chr1_rec_df,
                                     face_edge=('white', blue),
                                     mfrac_list=mfracs, flip=True)
    axrec = pI_by_targetmfrac_violin(axrec, chrX_rec_df,
                                     face_edge=('white', orange),
                                     mfrac_list=mfracs, flip=True)

    # ratios
    def extract_errorbar_data_e95ci(stats_dict):
        this_data = [stats_dict[m] for m in mfracs]
        these_means, _, _, these_centiles = zip(*this_data)

        lower_abs, upper_abs = zip(*these_centiles)
        lower_rel = these_means - np.asarray(lower_abs)
        upper_rel = np.asarray(upper_abs) - these_means
        yerr_matrix = np.vstack([lower_rel, upper_rel])

        return these_means, yerr_matrix

    def extract_errorbar_data_n95ci(stats_dict):
        this_data = [stats_dict[m] for m in mfracs]
        these_means, these_stderrs, these_ns, _ = zip(*this_data)

        yerr_list = [1.96 * sterr / n[0] for sterr, n in zip(these_stderrs, these_ns)]

        return these_means, yerr_list

    if errs == 'e':
        add_means, add_yerrs = extract_errorbar_data_e95ci(add_stats)
        rec_means, rec_yerrs = extract_errorbar_data_e95ci(rec_stats)
    else:
        add_means, add_yerrs = extract_errorbar_data_n95ci(add_stats)
        rec_means, rec_yerrs = extract_errorbar_data_n95ci(rec_stats)

    raxadd.errorbar(x=mfracs, y=add_means, yerr=add_yerrs, fmt='ko')
    raxrec.errorbar(x=mfracs, y=rec_means, yerr=rec_yerrs, fmt='ko', markerfacecolor='white')

    # lines
    # chrX sexbias initial introg lines
    initial_levels = 0.05 * np.array([0.5, 0.75, 1.0])
    levels = zip(initial_levels, mfracs)
    points = [((x - 0.2, x + 0.2), (i, i)) for i, x in levels]
    endpoints = list(sum(points, ()))
    axadd.plot(*endpoints, linestyle='--', color=orange, zorder=1, alpha=0.75)
    axrec.plot(*endpoints, linestyle='--', color=orange, zorder=1, alpha=0.75)
    # aut initial introg linesx
    initial_levels = 0.05 * np.array([1, 1, 1])
    levels = zip(initial_levels, mfracs)
    points = [((x - 0.2, x + 0.2), (i, i)) for i, x in levels]
    endpoints = list(sum(points, ()))
    axadd.plot(*endpoints, linestyle='--', color=blue, zorder=1, alpha=0.75)
    axrec.plot(*endpoints, linestyle='--', color=blue, zorder=1, alpha=0.75)
    # theoretic ratio lines
    initial_levels = [2 / (2 - m) for m in mfracs]
    levels = zip(initial_levels, mfracs)
    points = [((x - 0.2, x + 0.2), (i, i)) for i, x in levels]
    endpoints = list(sum(points, ()))
    raxadd.plot(*endpoints, linestyle='solid', color='y', zorder=1, alpha=0.75)
    raxrec.plot(*endpoints, linestyle='solid', color='y', zorder=1, alpha=0.75)

    # labels
    raxadd.set_ylabel('Aut:ChrX')
    raxadd.set_yticks([1, 2, 3])
    axadd.set_ylabel('Archaic coverage (per bp)')
    raxadd.set_title('Additive variants ($h=0.5$)', fontsize=10)
    raxrec.set_title('Recessive variants ($h=0.0$)', fontsize=10)

    # flip male fraction to female fraction
    raxadd.invert_xaxis()
    raxrec.invert_xaxis()

    # legend
    x_line = mlines.Line2D([0], [0], color=orange, label='Chromosome X')
    a_line = mlines.Line2D([0], [0], color=blue, label='Autosome')
    axadd.legend(handles=[a_line, x_line])

    # panel labels
    raxadd.text(0.04, 0.97 - 0.03, 'A', transform=raxadd.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='left')
    raxrec.text(0.04, 0.97 - 0.03, 'B', transform=raxrec.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='left')
    axadd.text(0.04, 0.97, 'C', transform=axadd.transAxes,
               fontsize=12, fontweight='bold', va='top', ha='left')
    axrec.text(0.04, 0.97, 'D', transform=axrec.transAxes,
               fontsize=12, fontweight='bold', va='top', ha='left')

    plt.tight_layout()
    pass


plot_sexbias_violins_with_ratio(additive_ratio_stats, recessive_ratio_stats)


# %%  Save figures

save_figures = False
if save_figures:
    plt.savefig('fig3_violins.svg',
                format='svg', dpi=600, bbox_inches='tight')
    plt.savefig('fig3_violins.png',
                format='png', dpi=600, bbox_inches='tight')
