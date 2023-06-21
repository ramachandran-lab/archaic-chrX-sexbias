#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 1.

    Summary of published empirical estimates of archaic coverage in global groups show greater coverage rates on autosomes than chromosome X.
        Each point indicates the per-bp amount of archaic coverage on autosomes and chromosome X in a global human group; see S1 Table for numerical values.
        Colors indicate continental group; see Table 1 for key.
        Shapes indicate estimation method; see S2 Table for details.
        Light grey clines show autosomal to chromosome X coverage ratios of 1 (i.e., equal per bp coverage), 2 (i.e., autosomal archaic coverage is double that of chromosome X), 4, 8, and 16.
        Dark grey cline shows median autosomal to chromosome X coverage ratio of 7.91.
        Nearly all points are well below the 1:1 line, indicating greater per-bp coverage on autosomes than chromosome X.
"""

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib import lines as mlines

# %% Import data

# Datafiles
# sgdp_info_file2 = 'fig1_data/10_24_2014_SGDP_metainformation_update.txt'
sgdp_info_file = 'fig1_data/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt'
sank_file = 'fig1_data/Sank16long.txt'
table1_file = 'fig1_data/rev2/table1_missing.txt'
skovSGDP_file = 'fig1_data/pgen.1007641.s012.txt'
hmmix_file = 'fig1_data/rev2/table1_hmmix.txt'

# Dataframes
oldtab1df = pd.read_csv(table1_file, sep=' ', index_col=False)
sankdf = pd.read_csv(sank_file, sep='\t', index_col=False, usecols=range(8),
                     names=['aut', 'aut_stderr', 'chrX', 'chrX_stderr',
                            'data_source', 'archaic', 'superpop', 'pop'])
sgdpdf = pd.read_csv(sgdp_info_file, sep='\t', index_col=False)
# sgdpdf2 = pd.read_csv(sgdp_info_file2, sep='\t', index_col=False)
hmmixdf = pd.read_csv(hmmix_file, sep=' ')


# %% Assemble data into one df

#  Assign continental codes across data sources
sp_dict = dict(zip(sgdpdf.Population_ID, sgdpdf.Region))
sp2sp_dict = dict(zip(['America', 'CentralAsiaSiberia', 'EastAsia', 'NA',
                       'Oceania', 'SouthAsia', 'WestEurasia'],
                      ['AMR', 'CAS', 'EAS', 'NA', 'OCE', 'SAS', 'EUR']))


def xlate_sp(region):
    return sp2sp_dict[region]


def find_sgdp_region(pop):
    if pop not in sp_dict.keys():
        if pop == 'Papuans':
            return 'Oceania'
        else:
            print(f'Can\'t find region for {pop}.')
        return 'NA'
    else:
        return sp_dict[pop]


def process_SGDP_Skov_data(datafile):
    all_data = pd.read_csv(datafile, sep='\t', header=0,
                           usecols=[0, 1, 4, 6, 8, 9, 14],
                           names=['names', 'chr', 'start', 'end', 'length',
                                  'num_snps', 'pop', 'country', 'region',
                                  'prob', 'snpA', 'snpD', 'snpV', 'outgroup',
                                  'method'])
    df = all_data[all_data['method'] == 'HMM']
    df = df[df['prob'] > 0.8]
    df = df.drop('prob', axis='columns')

    populations = df['pop'].unique()
    chromosomes = df['chr'].unique()
    samples = df['names'].unique()

    b = df.groupby(["pop", "chr"]).sum()
    c = df.groupby(["names", "chr"]).sum()

    # check that all pop'ns have all chromosomes
    assert all([p in b.index for p in zip(populations, chromosomes)])
    # check that all samples have all chromosomes
    assert all([p in c.index for p in zip(samples, chromosomes)])

    # There are 6230 entries in c, yet number of unique samples*23 = 6233
    [len(sdf["names"].unique()) for sdf in [df[df["chr"] == c] for c in chromosomes]]
    # three samples are missing chrX data.  which ones?
    aut_samps, x_samps = [sdf["names"].unique() for sdf in [df[df["chr"] == c]
                                                            for c in ['22', 'X']]]
    assert set(samples) == set(aut_samps)
    no_x_samps = list(set(aut_samps) - set(x_samps))
    # Remove the samples without chrX data from df
    df = df.drop(df[df["names"].isin(no_x_samps)].index)
    assert set(df["names"]) == set(x_samps)

    # Remove Papuans
    df = df[df['pop'] != 'Papuans']  # len 113444

    #  Make aut:chrX ratio
    # Collapse autosomes into sums
    a_sums = df[df["chr"] != "X"].groupby(['pop', 'names']).sum()
    x_sums = df[df["chr"] == "X"].groupby(['pop', 'names']).sum()
    a_sums.rename(columns={"length": "aut_tot"}, inplace=True)
    x_sums.rename(columns={"length": "chrX_tot"}, inplace=True)
    sums = pd.concat([a_sums, x_sums], axis=1)
    # Get lengths from Hg19 (aka GRCh37) to make per-bp numbers
    aut_len = 2981033286
    x_len = 155270560
    sums['aut'] = sums['aut_tot'] / aut_len
    sums['chrX'] = sums['chrX_tot'] / x_len

    df = sums.reset_index()
    # percentages to match other data
    df['aut'] *= 100
    df['chrX'] *= 100

    means = df.groupby('pop').mean()
    ses = df.groupby('pop').std()

    df = pd.merge(means, ses, how='outer', on=means.index,
                  suffixes=('', '_stderr'))
    df = df.drop(['aut_tot', 'chrX_tot', 'aut_tot_stderr', 'chrX_tot_stderr'],
                 axis='columns')
    df.rename(columns={"aut_mean": "aut", 'chrX_mean': 'chrX', 'key_0': 'pop'},
              inplace=True)
    df['data_source'] = 'SkovSGDP'
    return df


skovSGDPdf = process_SGDP_Skov_data(skovSGDP_file)

# %%
# Calculate ratios
sankdf['region'] = sankdf['pop'].apply(find_sgdp_region)
sankdf['superpop'] = sankdf['region'].apply(xlate_sp)
sank_ratios = sankdf["aut"] / sankdf["chrX"]
sankdf.insert(4, 'Aut_ChrX_ratio', sank_ratios)

skovSGDPdf['region'] = skovSGDPdf['pop'].apply(find_sgdp_region)
skovSGDPdf['superpop'] = skovSGDPdf['region'].apply(xlate_sp)
skovSGDP_ratios = skovSGDPdf['aut'] / skovSGDPdf['chrX']
skovSGDPdf.insert(4, 'Aut_ChrX_ratio', skovSGDP_ratios)

hmmixdf['data_source'] = 'hmmix'
hmmixdf['Aut_ChrX_ratio'] = hmmixdf['aut_x_ratio']

# Merge dataframes together
tab1df = pd.merge(oldtab1df, sankdf, how='outer')
tab1df = pd.merge(tab1df, skovSGDPdf, how='outer')
tab1df = pd.merge(tab1df, hmmixdf, how='outer')

superpop_colors = {'EAS': 'green', 'EUR': 'blue', 'AFR': 'DarkOrange',
                   'SAS': 'DarkMagenta', 'AMR': 'FireBrick',
                   'CAS': 'GoldenRod', 'OCE': 'DarkTurquoise'}
tab1df['spcol'] = [superpop_colors[sp] for sp in tab1df['superpop']]

tab1df['aut'] *= 0.01
tab1df['chrX'] *= 0.01


# %% Plotting

# params for figure
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

data_list = ['ARGWeaver-D', 'Dical-admix', 'Sank14', 'Sank16Neand',
             'Sank16Denis', 'hmmix', 'SkovSGDP']
marker_list = ['p', '*', 'o', '.', 'x', 'X', '+']
data2m_dict = dict(zip(data_list, marker_list))

fig, ax = plt.subplots(1, figsize=(7.5, 5))

# clines
tallest = max(tab1df['chrX'])
widest = max(tab1df['aut'])
allmedian = tab1df['Aut_ChrX_ratio'].median()
ax.plot([0, widest], [0, (1 / allmedian) * widest], ls='solid', color='grey', zorder=2)

ax.annotate("Aut:chrX::1:1", xy=(0.0138, 0.013))
ax.plot([0, tallest], [0, tallest], ls='solid', c='lightgrey', zorder=1)

ax.annotate("2:1", xy=(0.0272, 0.013))
ax.plot([0, 2 * tallest], [0, tallest], ls='solid', c='lightgrey', zorder=1)

ax.annotate("4:1", xy=(0.0279, 0.0074))
ax.plot([0, widest], [0, 0.25 * widest], ls='solid', c='lightgrey', zorder=1)

# ax.annotate("8:1", xy=(0.038, 0.0041))
# ax.plot([0, widest], [0, 0.125 * widest], ls='solid', c='lightgrey', zorder=1)

ax.annotate("16:1", xy=(0.02737, 0.001))
ax.plot([0, widest], [0, (1 / 16) * widest], ls='solid', c='lightgrey', zorder=1)

# data (scatter, no errorbars)
for d in data_list:
    mask = tab1df['data_source'] == d
    pltdf = tab1df[mask]
    if d == 'Sank16Neand':
        ax.scatter(pltdf['aut'], pltdf['chrX'], c=pltdf['spcol'],
                   marker=data2m_dict[d], label=d, alpha=1, s=5, zorder=3)
    else:
        ax.scatter(pltdf['aut'], pltdf['chrX'], c=pltdf['spcol'],
                   marker=data2m_dict[d], label=d, alpha=1, zorder=3)

ax.set_xlabel('Autosomal archaic coverage (per bp)')
ax.set_ylabel('X Chromosomal archaic coverage (per bp)')


# legends
# help with legend positioning https://stackoverflow.com/a/48405821/18621926
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
afr_patch = mpatches.Patch(edgecolor=superpop_colors['AFR'],
                           facecolor=superpop_colors['AFR'],
                           label='AFR')
oce_patch = mpatches.Patch(edgecolor=superpop_colors['OCE'],
                           facecolor=superpop_colors['OCE'],
                           label='OCE')
cas_patch = mpatches.Patch(edgecolor=superpop_colors['CAS'],
                           facecolor=superpop_colors['CAS'],
                           label='CAS')

sp_handles = [afr_patch, amr_patch, eas_patch, eur_patch, cas_patch, sas_patch, oce_patch]
sp_legend = plt.legend(title="Continental grouping", handles=sp_handles, frameon=True,
                       loc='lower left', bbox_to_anchor=(1, 0.5))
fig.add_artist(sp_legend)

marker_handles = [mlines.Line2D([], [], color='k', marker=m, linestyle='None',
                                label=d)
                  for m, d in zip(marker_list, data_list)]
ds_legend = plt.legend(title="Estimation method", handles=marker_handles,
                       frameon=True, loc='upper left', bbox_to_anchor=(1, 0.5))
fig.add_artist(ds_legend)


# Zoom axis inset
dax = ax.inset_axes([0.01, .65, 0.3, 0.3])
dax.plot([0, widest], [0, (1 / allmedian) * widest], ls='solid', color='grey', zorder=1)
dax.plot([0, tallest], [0, tallest], ls='solid', c='lightgrey', zorder=1)
dax.plot([0, 2 * tallest], [0, tallest], ls='solid', c='lightgrey', zorder=1)
dax.plot([0, widest], [0, 0.25 * widest], ls='solid', c='lightgrey', zorder=1)
# dax.plot([0, widest], [0, 0.125 * widest], ls='solid', c='lightgrey', zorder=1)
dax.plot([0, widest], [0, (1 / 16) * widest], ls='solid', c='lightgrey', zorder=1)
for d in data_list:
    mask = tab1df['data_source'] == d
    pltdf = tab1df[mask]
    if d == 'Sank16Neand':
        dax.scatter(pltdf['aut'], pltdf['chrX'], c=pltdf['spcol'],
                    marker=data2m_dict[d], label=d, alpha=1, s=5)
    else:
        dax.scatter(pltdf['aut'], pltdf['chrX'], c=pltdf['spcol'],
                    marker=data2m_dict[d], label=d, alpha=1)

dax.set_xlim([-0.0001, 0.0015])
dax.set_ylim([-0.0001, 0.0008])
dax.set_xticks([])
dax.set_yticks([])
# ax.indicate_inset_zoom(dax)


# %%  Save figures

save_figures = False
if save_figures:
    plt.savefig('fig1_rev2_table1.svg',
                format='svg', dpi=600, bbox_inches='tight')
    plt.savefig('fig1_rev2_table1.png',
                format='png', dpi=600, bbox_inches='tight')
