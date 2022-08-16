#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 5.

    The rapid timecourse of archaic coverage purging on the autosomes and chromosome X.

        Carets on the x-axis indicate the generation at which mean coverage has fallen halfway to its final value from its value in the generation immediately following introgression.
        Carets on the y-axis indicate the initial introgression fraction.
        Error bars are bootstrapped empirical 95% confidence intervals around means of all samples across all simulation replicates.
        See S5 Fig for a visualization the per-haplotype variance in archaic coverage over time.
    A
        Autosomal simulations.
    B
        Chromosome X simulations.
"""
import extract_pI_data as pI
import extract_pI_data_sexbias as sexpI

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as colors
from matplotlib.collections import PatchCollection
import numpy as np


# %% Import data

path_to_sex = 'fig5_data/'
sex_df = sexpI.compile_df_from_pI_files(path_to_sex, include_timecourse=True)

sex1_add_df = sexpI.subset_pI(sex_df, c='1', d='0.5SLiM', s=1, n=0.01, m=[0, 0.5, 1])
sex1_rec_df = sexpI.subset_pI(sex_df, c='1', d='0.0', s=1, n=0.01, m=[0, 0.5, 1])
sexX_add_df = sexpI.subset_pI(sex_df, c='X', d='0.5SLiM', s=1, n=0.01, m=[0, 0.5, 1])
sexX_rec_df = sexpI.subset_pI(sex_df, c='X', d='0.0', s=1, n=0.01, m=[0, 0.5, 1])


# %% Functions needed for plotting

def find_time_to_halfway(df, mf, from_f1=True):
    mfdf = sexpI.subset_pI(df, m=mf)
    timepoint_means = mfdf.groupby('timepoint')['pI'].mean()
    final_mean = timepoint_means[0]

    if from_f1:
        initial = timepoint_means.iloc[-1]
    else:
        initial = 0.05
    halfway = 0.5 * (initial - final_mean) + final_mean
    timepoint_diff_from_half = timepoint_means - halfway

    further_than_half = timepoint_diff_from_half[timepoint_diff_from_half.lt(0)]
    first_gen_further = further_than_half.tail(1).index.item()

    # simulations are run a double scale.  Change back to nominal generations.
    first_gen_simscale = 2 * 7500
    first_gen_further_simscale = 2 * first_gen_further

    time_to_half = first_gen_simscale - first_gen_further_simscale
    return time_to_half


def mean_pI_by_timepoint_errorbar(ax, input_df, label=None, style='o',
                                  color=None, fill=None, size=10, B=20,
                                  verbose=True):
    """Construct 95% confidence intervals reflecting grand samples.
    """
    def bs_mean(data, B=2000, thresh=5):
        n = len(data)
        mean = np.mean(data)
        draws = np.random.choice(data, size=(B, n), replace=True)
        bs_means = np.mean(draws, axis=1)
        ci_limits = np.percentile(bs_means, (thresh / 2, 100 - thresh / 2))
        yerr_lwr = mean - ci_limits[0]
        yerr_upr = ci_limits[1] - mean
        return mean, ci_limits, (yerr_lwr, yerr_upr)

    timepoints = np.unique(input_df['timepoint'])

    mean_list = []
    yerr_pair_list = []
    for t in timepoints:
        this_df = pI.subset_pI(input_df, t=t)
        if verbose:
            print(f"timepoint {t}/7500: n=={len(this_df['pI'])}, B=={B} ")
        timepoint_mean, this_e95CI, these_yerr_limits = bs_mean(this_df['pI'], B=B)
        mean_list.append(timepoint_mean)
        yerr_pair_list.append(these_yerr_limits)
    yerr_lwr, yerr_upr = zip(*yerr_pair_list)
    yerr_matrix = np.vstack([yerr_lwr, yerr_upr])

    # create "false"/"broken" x axis to accommodate timepoint 0 with false0
    broken_0 = timepoints[1] - 5  # x-axis location at which to place "0"
    broken_x_axis = [broken_0] + list(timepoints[1:])

    # call plot
    ax.errorbar(broken_x_axis, mean_list, yerr=yerr_matrix, label=label,
                fmt=style, mfc=fill, color=color, markersize=size, linewidth=size / 4)

    # label x-axis each 5 generations in reverse time
    ax.set_xlim([timepoints[-1] + 1, broken_0 - 1])
    timeticks = [int(i) for i in np.arange(start=timepoints[-1] + 1,
                                           stop=timepoints[1], step=-5)]
    ax.set_xticks(timeticks + [broken_0])
    tt_labels = [2 * t for t in timeticks] + [0]
    ax.set_xticklabels(tt_labels)

    return ax


# Make two-color boxes for legend
    # https://stackoverflow.com/a/67870930/18621926
# define an object that will be used by the legend
class MulticolorPatch(object):
    def __init__(self, colors):
        self.colors = colors


# define a handler for the MulticolorPatch object
class MulticolorPatchHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width, height = handlebox.width, handlebox.height
        patches = []
        for i, c in enumerate(orig_handle.colors):
            patches.append(plt.Rectangle([width / len(orig_handle.colors) * i - handlebox.xdescent,
                                          -handlebox.ydescent],
                           width / len(orig_handle.colors),
                           height,
                           facecolor=c,
                           edgecolor='none'))
        patch = PatchCollection(patches, match_original=True)
        handlebox.add_artist(patch)
        return patch


# %% Plot

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8


def plot_dom_timecourse_mean_by_sexbias_errorbar_2pack(x_add_df, x_rec_df,
                                                       a_add_df, a_rec_df,
                                                       size=10, verbose=False):
    fig, (axa, axx) = plt.subplots(ncols=2, sharey=True, sharex=True, figsize=(5.2, 3))

    # Legend
    cmap = colors.get_cmap(name='Set2')
    color_dict = {('X', 0): cmap(3), ('X', 0.5): cmap(1), ('X', 1): cmap(5),
                  ('1', 0): cmap(4), ('1', 0.5): cmap(0), ('1', 1): cmap(2)}
    male_colors = MulticolorPatch([color_dict[('1', 1)], color_dict[('X', 1)]])
    even_colors = MulticolorPatch([color_dict[('1', 0.5)], color_dict[('X', 0.5)]])
    female_colors = MulticolorPatch([color_dict[('1', 0)], color_dict[('X', 0)]])

    add_shape = mlines.Line2D([], [], color='k', marker='o', linestyle='None',
                              label='Additive variants ($h=0.5$)')
    rec_shape = mlines.Line2D([], [], color='k', marker='o', mfc='none',
                              linestyle='None', label='Recessive variants ($h=0.0$)')
    l = mlines.Line2D([0], [0], color="w")

    legend = fig.legend(handles=[add_shape,
                                 rec_shape,
                                 l,
                                 male_colors,
                                 even_colors,
                                 female_colors],
                        labels=['Additive variants ($h=0.5$)',
                                'Recessive variants ($h=0.0$)',
                                '',
                                'male introgressors ($p=0$)',
                                'no sex-bias ($p=0.5$)',
                                'female introgressors ($p=1$)'],
                        handler_map={MulticolorPatch: MulticolorPatchHandler()},
                        frameon=True, fontsize=8, ncol=2)

    # plot points
    def one_dom(input_ax, input_df, chromosome, marker='o', fill=None, verbose=False):
        mfrac_levels = input_df["target_male_frac"].unique()
        mfrac_levels.sort()
        for mf in mfrac_levels:
            mdf = sexpI.subset_pI(input_df, m=mf)
            col = color_dict[(chromosome, mf)]

            input_ax = mean_pI_by_timepoint_errorbar(input_ax, mdf,
                                                     label='p = ' + str(1 - mf),
                                                     color=col, style=marker,
                                                     fill=fill, size=size,
                                                     B=2000, verbose=verbose)
            # x-axis carets for 50% timepoint
            half_x = 7500 - 0.5 * find_time_to_halfway(input_df, mf) - 0.15 + 0.3 * mf
            input_ax.plot(half_x, 0, marker=10, color=col, mfc=fill, clip_on=False)

            # y-axis carets for initial introg line
            level = 0.05
            if chromosome == 'X':
                level *= (1 - 0.5 * mf)
            input_ax.plot(7500, level, marker=9, color=col, mfc=fill, clip_on=False)
        return input_ax

    if verbose:
        print('---Generating errorbars for chrX additive---')
    axx = one_dom(axx, x_add_df, 'X', marker='o', verbose=verbose)
    if verbose:
        print('---Generating errorbars for chrX recessive---')
    axx = one_dom(axx, x_rec_df, 'X', marker='o', fill='none', verbose=verbose)
    if verbose:
        print('---Generating errorbars for aut additive---')
    axa = one_dom(axa, a_add_df, '1', marker='o', verbose=verbose)
    if verbose:
        print('---Generating errorbars for aut recessive---')
    axa = one_dom(axa, a_rec_df, '1', marker='o', fill='none', verbose=verbose)

    # labels and title
    axa.set_ylim([0, 0.065])
    axa.set_ylabel('Archaic coverage (per bp)')
    axx.set_xlabel('Generations ago')
    axa.set_xlabel('Generations ago')

    axa.set_title('Autosome')
    axx.set_title('Chromosome X')

    for ax, label in zip(fig.get_axes(), ['A', 'B', 'C', 'D', 'E']):
        ax.text(0.1, 0.97, label, transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='right')

    plt.tight_layout()
    pass


#  Generating errorbars for all points can be slow.
#  Set verbose=True for progress indication.
plot_dom_timecourse_mean_by_sexbias_errorbar_2pack(sexX_add_df, sexX_rec_df,
                                                   sex1_add_df, sex1_rec_df,
                                                   size=4, verbose=False)

# %% Save

# plt.savefig('fig5_timecourse.svg',
#             format='svg', dpi=600, bbox_inches='tight')
# plt.savefig('fig5_timecourse.png',
#             format='png', dpi=600, bbox_inches='tight')
