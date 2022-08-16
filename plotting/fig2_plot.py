#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 2.

    Distribution of archaic coverage between autosomes and chromosome X amongst 1kG samples.
        Coverage estimated by Skov et al. [18] from data from The 1000 Genomes Project Consortium [19].
        See Table 1 for key to abbreviations.
    A
        Mean aut:chrX archaic coverage ratio among continental groups.
        Error bars show bootstrapped 95% confidence intervals around the means.
        Distributions of archaic coverage on autosomes
    B
        and chromosome X
    C
        across individuals, colored by continental group. Carets indicate means.
        Summaries of coverage ratios
    D,
        as well as coverage distributions on the two chromosome types
    E
        among 1kG sample groups.
        Populations are sorted by ascending archaic coverage amounts on chromosome X.
"""

import pandas as pd
import numpy as np
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


# %% Assemble data (Skov 1kG calls)

datafile_1kG = 'fig2_data/LauritsSkovIntrogressionMaps/Archaicsegment_1000genomes.txt.gz'

all_data_1kG = pd.read_csv(datafile_1kG, sep='\t', header=0,
                           names=['names', 'chr', 'start', 'end',
                                  'length', 'num_snps', 'pop', 'prob',
                                  'snpA', 'snpD', 'snpV', 'snps'],
                           usecols=[0, 1, 4, 6, 7])

df = all_data_1kG[all_data_1kG['prob'] > 0.8]

df = df.drop('prob', axis='columns')
populations = df['pop'].unique()
chromosomes = df['chr'].unique()
samples = df['names'].unique()
# our populations are:
#    ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL',
#     'KHV', 'BEB', 'STU', 'ITU', 'CEU', 'CHB', 'JPT', 'MXL', 'TSI',
#     'GIH']

b = df.groupby(["pop", "chr"]).sum()
c = df.groupby(["names", "chr"]).sum()

# check that all pop'ns have all chromosomes
assert all([p in b.index for p in zip(populations, chromosomes)])

# check that all samples have all chromosomes
assert all([p in c.index for p in zip(samples, chromosomes)])

# There are 42385 entries in c, yet number of unique samples*23 = 42389
[len(sdf["names"].unique()) for sdf in [df[df["chr"] == c] for c in chromosomes]]
# four samples are missing chrX data.  which ones?
aut_samps, x_samps = [sdf["names"].unique() for sdf in [df[df["chr"] == c] for
                                                        c in ['22', 'X']]]
assert all(aut_samps == samples)
no_x_samps = list(set(aut_samps) - set(x_samps))
# Remove the samples without chrX data from df
df = df.drop(df[df["names"].isin(no_x_samps)].index)
assert set(df["names"]) == set(x_samps)

# #  Make aut:chrX ratio
# Collapse autosomes into sums
a_sums = df[df["chr"] != "X"].groupby(['pop', 'names']).sum()
x_sums = df[df["chr"] == "X"].groupby(['pop', 'names']).sum()

a_sums.rename(columns={"length": "aut_len"}, inplace=True)
x_sums.rename(columns={"length": "x_len"}, inplace=True)
sums = pd.concat([a_sums, x_sums], axis=1)

# Get lengths from Hg19 (aka GRCh37) to make per-bp numbers
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
aut_len = 2981033286
x_len = 155270560

sums['aut_pbp'] = sums['aut_len'] / aut_len
sums['x_pbp'] = sums['x_len'] / x_len

sums['aut_x_ratio'] = sums['aut_pbp'] / sums['x_pbp']


# # make a df
dfp = sums.reset_index(level=['pop'])
superpops = dfp["pop"].apply(find_superpop)
dfp.insert(0, 'superpop', superpops)
# remove samples with outlying aut:chrX ratios
midf = dfp[dfp['aut_x_ratio'] < 200]

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


# %% **Multipanel figure
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10


fig = plt.figure(figsize=(7.5, 7.5))
gs = fig.add_gridspec(nrows=4, ncols=2)

axa = fig.add_subplot(gs[0:2, 0])
axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[1, 1], sharex=axb)
axd = fig.add_subplot(gs[2, :])
axe = fig.add_subplot(gs[3, :])


# # Panels A, B, C
plot_ratio_eCI(midf, by='superpop', existing_axis=axa)
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
    caret_y = 80
    axb.plot(aut_means[-1], caret_y, marker=7, clip_on=False, ls='none',
             markerfacecolor=superpop_colors[sp], markeredgecolor=aec)
    axc.plot(x_means[-1], caret_y, marker=7, clip_on=False, ls='none',
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
axe.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
axe.tick_params(axis='x', pad=12)

# Panel D: Ratio errorbar
axd.set_ylabel("Aut:ChrX")
axd = plot_ratio_eCI(midf, 'pop', sort_order=popns, existing_axis=axd)
axd.tick_params(labelbottom=False)
axd.set_xlim(axe.get_xlim())

# # Legend
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
for ax, label in zip(fig.get_axes(), ['A', 'B', 'C', 'D', 'E']):
    ax.text(0.04, 0.97, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')

plt.savefig('fig2_skov1kG_1labels.png',
            format='png', dpi=600, bbox_inches='tight')
# plt.savefig('fig2_skov1kG.svg',
#             format='svg', dpi=600, bbox_inches='tight')
