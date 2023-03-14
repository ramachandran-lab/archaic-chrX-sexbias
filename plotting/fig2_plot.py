#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 2.

    Distribution of archaic coverage between autosomes and chromosome X amongst 1kG samples.
        Coverage estimated by method of Skov et al. [18] from data from The 1000 Genomes Project Consortium [19].
        See Table 1 for key to abbreviations.
    A
        Mean aut:chrX archaic coverage ratio among continental groups.
        Error bars show bootstrapped 95% confidence intervals around the means.
    B
        Distributions of archaic coverage on autosomes and
    C
        chromosome X across individuals, colored by continental group.
        Carets indicate means.
    D
        Summaries of coverage ratios, as well as
    E
        coverage distributions on the two chromosome types among 1kG sample groups.
        Populations are sorted by ascending archaic coverage amounts on chromosome X.
    F
        Rank orderings of 1kG sample groups by mean chrX coverage, mean autosomal coverage, and aut:chrX ratio, sorted from least to greatest.
"""

import pandas as pd
import numpy as np
from interlap import InterLap
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.cm as colors


def find_superpop(pop, return_color=False):
    kG_popns = ['CHB',  'JPT',  'CHS',  'CDX',  'KHV',  'CEU',  'TSI',  'FIN',
                'GBR',  'IBS',  'YRI',  'LWK',  'MAG',  'MSL',  'ESN',  'ASW',
                'ACB',  'MXL',  'PUR',  'CLM',  'PEL',  'GIH',  'PJL',  'BEB',
                'STU',  'ITU']
    kG_superpopns = ['EAS',  'EAS',  'EAS',  'EAS',  'EAS',  'EUR',  'EUR',
                     'EUR',  'EUR',  'EUR',  'AFR',  'AFR',  'AFR',  'AFR',
                     'AFR',  'AFR',  'AFR',  'AMR',  'AMR',  'AMR',  'AMR',
                     'SAS',  'SAS',  'SAS',  'SAS',  'SAS']

    idx = kG_popns.index(pop)
    superpop = kG_superpopns[idx]

    if return_color:
        c = superpop_colors[superpop]
        return c
    else:
        return superpop


# %% Grab data from hmmix
hmmix_file = f'fig2_data/noPARnoselection_output.csv.gz'

df = pd.read_csv(hmmix_file)
df['names'] = df['name']
# remove samples
midf = df[df['aut_pbp'] < 0.1]
midf = midf[midf['aut_x_ratio'] < 50]

# %%  create file for table 1

tab1_means = midf.groupby('pop').mean()[['aut_pbp', 'x_pbp', 'aut_x_ratio']]
tab1_means['aut'] = tab1_means['aut_pbp']
tab1_means['aut'] *= 100
tab1_means['aut'] = tab1_means['aut'].apply(lambda x: round(x, 2))

tab1_means['chrX'] = tab1_means['x_pbp']
tab1_means['chrX'] *= 100
tab1_means['chrX'] = tab1_means['chrX'].apply(lambda x: round(x, 2))

tab1_means['aut_x_ratio'] = tab1_means['aut_x_ratio'].apply(lambda x: round(x, 1))
tab1_means.drop(columns=['aut_pbp', 'x_pbp'], inplace=True)

tab1_sterr = midf.groupby('pop').std()[['aut_pbp', 'x_pbp']]
tab1_sterr *= 100
tab1_sterr = tab1_sterr.apply(lambda x: round(x, 2))

tab1 = tab1_means.join(tab1_sterr)
tab1['aut_stderr'] = tab1['aut_pbp']
tab1['chrX_stderr'] = tab1['x_pbp']
tab1.drop(columns=['aut_pbp', 'x_pbp'], inplace=True)
tab1.sort_values('aut_x_ratio', inplace=True, ascending=False)

tab1.to_csv('table1_hmmix.txt', sep=' ')

# %% Functions needed for plotting

# colors
cmap = colors.get_cmap(name='tab20')
xec = cmap(2)
aec = cmap(0)
superpop_colors = {'EAS': 'green', 'EUR': 'blue', 'AFR': 'DarkOrange',
                   'SAS': 'DarkMagenta', 'AMR': 'FireBrick'}
sps = ['AMR', 'EAS', 'EUR', 'SAS']
spcolors = [superpop_colors[sp] for sp in sps]


# Change sort order of populations for panels D, E
def get_sort_order(by, df=midf):
    sdf = midf.sort_values(by='pop')
    means_df = sdf.groupby(['superpop', 'pop', 'names']).mean()
    # sort
    gb = means_df[by].groupby('pop')
    means = gb.mean()
    means = means.sort_values(ascending=True)
    sort_order = means.index
    return sort_order


popns = get_sort_order('x_pbp')


def plot_ratio_eCI(df, by, sort_order=None, existing_axis=None, thresh=5):

    def bs_mean(data, B=2000, thresh=5):
        n = len(data)
        mean = np.mean(data)
        draws = np.random.choice(data, size=(B, n), replace=True)
        bs_means = np.mean(draws, axis=1)
        ci_limits = np.percentile(bs_means, (thresh / 2, 100 - thresh / 2))
        yerr_lwr = mean - ci_limits[0]
        yerr_upr = ci_limits[1] - mean
        return mean, ci_limits, (yerr_lwr, yerr_upr)

    if sort_order is None:
        sort_order = np.unique(df[by])

    means = []
    yerr_pairs = []
    for j in sort_order:
        jdf = df[df[by] == j]
        this_mean, _, these_yerrs = bs_mean(jdf['aut_x_ratio'])
        means.append(this_mean)
        yerr_pairs.append(these_yerrs)
    lower_errs, upper_errs = zip(*yerr_pairs)

    if existing_axis is None:
        fig, existing_axis = plt.subplots(1)
        existing_axis.set_ylabel('Autosomal:ChrX archaic coverage ratio')

    if by == 'superpop':
        color_list = [superpop_colors[sp] for sp in sort_order]
    elif by == 'pop':
        color_list = [find_superpop(p, return_color=True) for p in sort_order]
    else:
        color_list = 'k'

    existing_axis.errorbar(sort_order, means, yerr=[lower_errs, upper_errs],
                           fmt='none', ecolor=color_list)
    existing_axis.scatter(sort_order, means, c=color_list)
    return existing_axis



# %% * Plot Fig 2

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10


fig = plt.figure(figsize=(7.5, 7.5))
gs = fig.add_gridspec(nrows=5, ncols=2)

axa = fig.add_subplot(gs[0:2, 0])
axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[1, 1], sharex=axb)
axd = fig.add_subplot(gs[2, :])
axe = fig.add_subplot(gs[3, :])
axf = fig.add_subplot(gs[4, :])


# # Panels A, B, C
# Panel A:
plot_ratio_eCI(midf, by='superpop', existing_axis=axa, sort_order=['EUR', 'AMR', 'EAS', 'SAS'])
axa.set_ylabel('Autosomal : Chromosome X coverage ratio')

# Panels B, C: stacked bar by superpop
num_bins = 75

rel_x_by_superpop = []
rel_a_by_superpop = []
aut_means = []
x_means = []
for sp in sps:
    spdf = midf[midf['superpop'] == sp]
    x_means.append(np.mean(spdf['x_pbp']))
    aut_means.append(np.mean(spdf['aut_pbp']))
    rel_x_by_superpop.append(np.asarray(spdf['x_pbp']))
    rel_a_by_superpop.append(np.asarray(spdf['aut_pbp']))

    # plot carets in loop for coloring
    axb_caret_y = 60
    axc_caret_y = 82.95
    axb.plot(aut_means[-1], axb_caret_y, marker=7, clip_on=False, ls='none',
             markerfacecolor=superpop_colors[sp], markeredgecolor=aec)
    axc.plot(x_means[-1], axc_caret_y, marker=7, clip_on=False, ls='none',
             markerfacecolor=superpop_colors[sp], markeredgecolor=xec)

axb.hist(rel_a_by_superpop, bins=num_bins, stacked=True,
         color=spcolors, label=sps)
axc.hist(rel_x_by_superpop, bins=num_bins, stacked=True,
         color=spcolors, label=sps)
axb.tick_params(labelbottom=False)
axc.set_xlabel('Archaic coverage per base pair')

# color axis boxes by aut/chrX
for spine in axc.spines.values():
    spine.set_edgecolor(xec)
for spine in axb.spines.values():
    spine.set_edgecolor(aec)

axc.tick_params(color=xec, labelcolor=xec, axis='y')
axb.tick_params(color=aec, labelcolor=aec)

for ticks in axc.xaxis.get_ticklines() + axc.yaxis.get_ticklines():
    ticks.set_color(xec)
for ticks in axb.xaxis.get_ticklines() + axb.yaxis.get_ticklines():
    ticks.set_color(aec)

axb.set_title('Autosomes', fontsize=10)
axc.set_title('Chromosome X', fontsize=10)

# # Panels D, E

# Panel E: Boxes
boxdata_x = [np.asarray(midf[midf['pop'] == p]['x_pbp']) for p in popns]
boxdata_a = [np.asarray(midf[midf['pop'] == p]['aut_pbp']) for p in popns]
# autosomes
pa_boxes = axe.boxplot(boxdata_a, patch_artist=True, positions=list(range(len(popns))))
# https://stackoverflow.com/a/20291461/18621926
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(pa_boxes[element], color=aec)
plt.setp(pa_boxes['fliers'], markeredgecolor=aec)
for patch, pop in zip(pa_boxes['boxes'], popns):
    patch.set_facecolor(find_superpop(pop, return_color=True))
# chromosome X
px_boxes = axe.boxplot(boxdata_x, patch_artist=True, positions=list(range(len(popns))))
for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
    plt.setp(px_boxes[element], color=xec)
plt.setp(px_boxes['fliers'], markeredgecolor=xec)
for patch, pop in zip(px_boxes['boxes'], popns):
    patch.set_facecolor(find_superpop(pop, return_color=True))
axe.set_xticks(list(range(len(popns))), labels=popns, fontsize=9)
axe.set_ylabel("Archaic coverage (per bp)")
axe.tick_params(top=True, labeltop=True, bottom=True, labelbottom=True)
axe.tick_params(axis='x', pad=12)

# Panel D: Ratio errorbar
axd.set_ylabel("Aut:ChrX")
axd = plot_ratio_eCI(midf, 'pop', sort_order=popns, existing_axis=axd)
axd.tick_params(labelbottom=False)
axd.set_xlim(axe.get_xlim())

# Panel F: Rankings
popns_x = get_sort_order('x_pbp')
popns_rat = get_sort_order('aut_x_ratio')
popns_aut = get_sort_order('aut_pbp')

rank_x = dict(zip(popns_x, range(len(popns_x))))
rank_aut = dict(zip(popns_aut, range(len(popns_aut))))
rank_rat = dict(zip(popns_rat, range(len(popns_rat))))

axf.set_ylim([-0.2, 2.2])
axf.set_yticks([0, 1, 2])
axf.set_yticklabels(['Aut:ChrX', 'aut\ncoverage', 'chrX\ncoverage'],
                    fontsize=9)
axf.set_xticks(list(range(len(popns_x))))
axf.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)
axf.set_xlim(axe.get_xlim())
axfb = axf.secondary_xaxis('bottom')
axfb.set_xticks([0, 18], labels=['least', 'greatest'])
for pop in popns_x:
    axf.plot([rank_x[pop], rank_aut[pop], rank_rat[pop]], [2, 1, 0],
             ls='solid', alpha=0.1,
             color=find_superpop(pop, return_color=True))
    marker_string = '$\mathtt{' + str(pop) + '}$'
    axf.plot([rank_x[pop], rank_aut[pop], rank_rat[pop]], [2, 1, 0],
             marker=marker_string, markersize=16, ls='None',
             color=find_superpop(pop, return_color=True))


# Legend
eas_patch = mpatches.Patch(edgecolor=superpop_colors['EAS'],
                           facecolor=superpop_colors['EAS'],
                           label='EAS')
eur_patch = mpatches.Patch(edgecolor=superpop_colors['EUR'],
                           facecolor=superpop_colors['EUR'],
                           label='EUR')
sas_patch = mpatches.Patch(edgecolor=superpop_colors['SAS'],
                           facecolor=superpop_colors['SAS'],
                           label='SAS')
amr_patch = mpatches.Patch(edgecolor=superpop_colors['AMR'],
                           facecolor=superpop_colors['AMR'],
                           label='AMR')
x_line = mlines.Line2D([0], [0], color=xec, label='Chromosome X')
a_line = mlines.Line2D([0], [0], color=aec, label='Autosome')

sp_handles = [amr_patch, eas_patch, eur_patch, sas_patch]
chr_handles = [a_line, x_line]
l = mlines.Line2D([0], [0], color="w")
legend_handles = sp_handles + [l] + chr_handles + [l]
axa.legend(handles=legend_handles, ncol=2, frameon=True)

plt.tight_layout()

# panel labels
for ax, label in zip(fig.get_axes(), ['A', 'B', 'C', 'D', 'E', 'F']):
    ax.text(-0.015, 1.05, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')

# %%  Save figures

save_figures = False
if save_figures:
    plt.savefig(f'fig2_hmmix.png', format='png', dpi=600, bbox_inches='tight')
    plt.savefig(f'fig2_hmmix.svg', format='svg', dpi=600, bbox_inches='tight')
