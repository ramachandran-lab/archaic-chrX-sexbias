#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting Fig 4.

    Archaic coverage tracts on chromosome X are longer than those on an autosome.
        All results shown reflect an additive model of dominance.
    A
        Tract length spectra from simulations of an archaic introgression scenario under autosomal or chromosome X inheritance; note log scale.
        Points indicate frequencies of archaic coverage tracts within length bins of 3000 bp.
        Simulated chromosomes shared a variable recombination rate landscape (see Methods for details).
        Solid lines indicate expected archaic coverage lengths under a simple, neutral model with constant recombination rate equal to the mean rate of the simulated recombination landscape.
    B
        Average recombination rate within each simulated archaic coverage tract, plotted against length of the tract.  Note log scale.
    C
        Proportion of archaic coverage tracts that are a given length or shorter; proportions are calculated within each inheritance type.
        Arrow indicates 95th percentile of tract lengths, illustrated in panel D alongside sex-biased scenarios.
    D
        Distribution of 95th percentiles of archaic coverage tract lengths found on either an autosome or chromosome X, across degrees of introgressor sex-bias.
        Each box reflects the length distributions from ten simulation replicates.
"""

from matplotlib import pyplot as plt
import matplotlib.patches as mpatch
import matplotlib.lines as mlines
import matplotlib.cm as colors

import numpy as np
import itertools


# %%  Go get the data
# Run prep_fig4.py to generate these pickle files.

import pickle

with open('fig4_tracts_A.pkl', "rb") as f:
    x_ax = pickle.load(f)
    a_counts, x_counts = pickle.load(f)
    exp_cdfA, exp_cdfX = pickle.load(f)

with open('fig4_tracts_B.pkl', "rb") as f:
    a_len_rec = pickle.load(f)
    x_len_rec = pickle.load(f)

with open('fig4_tracts_C.pkl', "rb") as f:
    a_props, x_props = pickle.load(f)

with open('fig4_tracts_D.pkl', "rb") as f:
    boxplot_cutoffs_add = pickle.load(f)


# %% Plot

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8

# colors
cmap = colors.get_cmap(name='tab20')
afc = cmap(1)
aec = cmap(0)
xfc = cmap(3)
xec = cmap(2)


fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(7.5, 5))

axa = axes[0][0]
axb = axes[1][0]
axc = axes[0][1]
axd = axes[1][1]

axa.get_shared_x_axes().join(axa, axb)
axa.tick_params(labelbottom=False)

# axis labels
axa.set_ylabel('Frequency')
axb.set_ylabel('Average recombination rate')
axc.set_ylabel('Cumulative proportion of tracts')
axd.set_ylabel('95%tile of coverage tract lengths (kbp)')

axb.set_xlabel('Length of archaic coverage tract (bp)')
axc.set_xlabel('Length of archaic coverage tract (bp)')


# panel labels
axa.text(0.98, 0.97, 'A', transform=axa.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='right')
axb.text(0.98, 0.97, 'B', transform=axb.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='right')
axc.text(0.06, 0.97, 'C', transform=axc.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='right')
axd.text(0.06, 0.97, 'D', transform=axd.transAxes,
         fontsize=12, fontweight='bold', va='top', ha='right')

plt.tight_layout(pad=0.5)


# Panel A
axa.plot(x_ax, a_counts[1:], marker='.', linestyle='None', label="Autosome")
axa.plot(x_ax, x_counts[1:], marker='.', linestyle='None', label="Chromosome X")
axa.set_ylim([1, 10e6])
axa.plot(x_ax, exp_cdfA, color=aec, alpha=0.8)
axa.plot(x_ax, exp_cdfX, color=xec, alpha=0.8)
axa.set_yscale('log')
axa.set_xscale('log')

#   legend
exp_line = mlines.Line2D([], [], color='k', marker=None, linestyle='solid',
                         alpha=0.8,
                         label='Neutral, constant recomb')
blue_blob = mpatch.Patch(color=aec, label='Autosome')
orange_blob = mpatch.Patch(color=xec, label='Chromosome X')
sim_dot = mlines.Line2D([], [], color='k', marker='.',
                        linestyle='None', label='Simulated tract lengths')

legend_info = [blue_blob, orange_blob, sim_dot, exp_line]
axa.legend(handles=legend_info, fontsize=8)

# Panel C
axc.plot(x_ax, a_props, marker='.', linestyle=':')
axc.plot(x_ax, x_props, marker='.', linestyle=':')
axc.set_xscale('log')
axc.annotate("", xy=(36000, 0.95), xytext=(20000, 0.95), arrowprops=dict(arrowstyle="->"))

# Panel B
a_len, a_rec = a_len_rec
x_len, x_rec = x_len_rec
a_len = np.concatenate(a_len).ravel()
x_len = np.concatenate(x_len).ravel()
a_rec = list(itertools.chain(*a_rec))
x_rec = list(itertools.chain(*x_rec))
axb.plot(a_len, a_rec, marker='.', lw=0, linestyle="", markersize=1,
         c=aec, rasterized=True)
axb.plot(x_len, x_rec, marker='.', lw=0, linestyle="", markersize=1,
         c=xec, rasterized=True)
axb.set_xscale('log')
axb.set_yscale('log')

# Panel D
bp_co_addA, bp_co_addX = zip(*boxplot_cutoffs_add)
kbp_co_addA = [box_data / 1000 for box_data in bp_co_addA]
kbp_co_addX = [box_data / 1000 for box_data in bp_co_addX]

cmfs = np.asarray([0, 0.5, 1])
scootch = 0.05
axd.boxplot(kbp_co_addA,
                positions=cmfs + scootch,
                widths=scootch * 1,
                patch_artist=True,
                boxprops=dict(facecolor=afc, color=aec),
                capprops=dict(color=aec),
                whiskerprops=dict(color=aec),
                flierprops=dict(color=afc, markeredgecolor=aec),
                medianprops=dict(color=aec)
            )

axd.boxplot(kbp_co_addX,
                positions=cmfs - scootch,
                widths=scootch * 1,
                patch_artist=True,
                boxprops=dict(facecolor=xfc, color=xec),
                capprops=dict(color=xec),
                whiskerprops=dict(color=xec),
                flierprops=dict(color=xfc, markeredgecolor=xec),
                medianprops=dict(color=xec)
            )
axd.set_xticks(cmfs)
axd.set_xticklabels(['female\nintrogressors\n($p = 1$)',
                     'no\nsexbias\n($p = 0.5$)',
                     'male\nintrogressors\n($p = 0$)'],
                    fontsize=8)
axd.invert_xaxis()


# %% Save

# plt.savefig('fig4_tracts.png',
#             format='png', dpi=600, bbox_inches='tight')
# plt.savefig('fig4_tracts.svg',
#             format='svg', dpi=600, bbox_inches='tight')
