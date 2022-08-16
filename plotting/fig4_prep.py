#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for preparing data for Fig 4.

See plot_fig4.py for plotting.
"""
from extract_pI_data_sexbias import get_properties_from_filename
from extract_pI_data_sexbias import subset_pI as subset_df

import glob
from os import path
import pandas as pd
import numpy as np
from interlap import InterLap

import pickle


def build_filename_df(path_to_files, no_ut=False):
    paths = glob.iglob(path_to_files + '*.npz')
    files = [path.basename(p) for p in paths]
    df = pd.DataFrame(columns=['chromosome', 'dominance', 's_scale', 'NsNr',
                               'target_male_frac', 'seed', 'filename'])
    for i, f in enumerate(files):
        chrom, dom, sscale, nsnr, mfrac, seed = get_properties_from_filename(f, '.tracts.npz', no_ut=no_ut)
        df.loc[i] = [chrom, dom, sscale, nsnr, mfrac, seed, path_to_files + f]
    return df


tract_files = 'fig4_data/'
df = build_filename_df(tract_files)
df_no_sexbias = subset_df(df, m=0.5, s=1)


# Save data needed to run `plot_fig4.py` as `.pkl` files.
save_fig4_data = True


# %% Panels A and B

# Parameters

region_len = 100000000  # simulating 100 megabase region
chr_info_file = '../generating_simulation_data/sim_seq_info.txt'
dominance = '0.5SLiM'  # additive variants
bin_width = 3000  # bin coverage tracts into bins of 3000 base pairs
m = 0.05  # 5% introgression pulse
num_haps = 1000  # 1000 haplotype samples from each simulation replicate
num_reps = 10  # 10 simulation replicates per scenario
# length bins for counting/summarizing coverage tract length spectra
longest_tract = 916181
bins = list(range(0, longest_tract + bin_width + 1, bin_width))
bin_midpoints = [b + (bin_width / 2) for b in bins]
x_ax = bin_midpoints[:-1]
num_bins = len(bins)


# Panel A: Calculate expected tract length spectra

g = 7500 * 2  # generations between introgression event and sampling
r_A = 9.764035057e-09  # average simulated autosomal per-bp recombination rate
#   in SLiM:  g = 7500, r = 0.5*(1-(1-2*9.764035057e-09)^2) = 1.952806999439005e-08
#   These are essentially the same.  99.99999938749734%
r_X = r_A * (2 / 3)  # scaled bulk chromosome X recombination rate
rate_A = g * r_A  # autosomal rate parameter for exponential tract breakage model
rate_X = g * r_X


def get_expected_num_blocks(rate, m, L, num_haps, num_reps):

    num_blocks_per_hap = L * rate
    num_blocks = num_blocks_per_hap * num_haps * num_reps
    num_neand_blocks = m * num_blocks
    return num_neand_blocks


expected_num_X = get_expected_num_blocks(rate_X, m=m, L=region_len,
                                         num_haps=num_haps, num_reps=num_reps)
expected_num_A = get_expected_num_blocks(rate_A, m=m, L=region_len,
                                         num_haps=num_haps, num_reps=num_reps)


def get_spectrum_expectation(num_tracts, bins, rate):
    def cdfexp(x, rate):
        return 1 - np.exp(-rate * x)

    def cdfdiff(right, left, rate):
        return cdfexp(right, rate) - cdfexp(left, rate)

    lbins = bins[:-1]
    rbins = bins[1:]

    cdfchunks = np.asarray([cdfdiff(r, l, rate) for r, l in zip(rbins, lbins)])
    expected_num_tracts = num_tracts * cdfchunks
    return expected_num_tracts


exp_cdfX = get_spectrum_expectation(expected_num_X, bins, rate_X)
exp_cdfA = get_spectrum_expectation(expected_num_A, bins, rate_A)


# Panel B: Collect recombination rate info used in simulations

def make_region_dfs(tract_map_fn):
    df = pd.read_csv(tract_map_fn, delim_whitespace=True,
                     names=('feature', 'data1', 'data2'))

    # separate into dfs and rename columns
    dfe = df[df['feature'] == 'exon']
    dfr = df[df['feature'] == 'recRate']
    dfe = dfe.rename({'feature': 'feature', 'data1': 'start', 'data2': 'end'},
                     axis='columns')
    dfr = dfr.rename({'feature': 'feature', 'data1': 'pos', 'data2': 'rate'},
                     axis='columns')

    # exon density
    dfe['length'] = dfe['end'] - dfe['start']
    exon_density = dfe['length'].sum() / region_len
    print(f"Exon density is {exon_density}.")

    # avg. rec rate
    pos_start = np.asarray(dfr['pos'][:-1])
    pos_start = np.insert(pos_start, 0, 0)
    dfr['length'] = dfr['pos'] - pos_start
    avg_rec_rate = np.average(dfr['rate'], weights=dfr['length'])
    print(f"Average recombination rate is {avg_rec_rate}.")

    # add start and end positions to recrates
    dfr['start'] = dfr['pos'] - dfr['length'] + 1
    dfr['end'] = dfr['pos']

    return dfe, dfr


def make_recrate_interlap(dfrr):
    inter = InterLap()
    dfrr = dfrr[['start', 'end', 'rate']]
    tuples = [tuple(x) for x in dfrr.to_records(index=False)]
    inter.update(tuples)
    return inter


_, df_recrates = make_region_dfs(chr_info_file)
recrate_interlap = make_recrate_interlap(df_recrates)


# Panels A&B: Count tracts in simulated data

def count_samplewise_tracts_with_recrates(df, chromosome, bins,
                                          recrate_interlap,
                                          dominance=dominance):
    def union_tracts(a):
        # modified from https://stackoverflow.com/a/15273749
        b = []
        for begin, end in sorted(a):
            if b and b[-1][1] >= begin - 1:
                b[-1][1] = max(b[-1][1], end)
            else:
                b.append([begin, end])
        return b

    def get_homwise_endpoints(endpoints_per_pairing):
        endpoints = sum(endpoints_per_pairing, [])  # flatten list of endpoint pairs
        tracts = union_tracts(endpoints)
        return tracts

    def import_samplewise_info(filepath):
        with np.load(filepath, allow_pickle=True) as all_spectra:
            # tract endpoints
            endpoints_pp = all_spectra['homwise_info'][:, 1]
            samplewise_tract_endpoints = [get_homwise_endpoints(pp) for
                                          pp in endpoints_pp]
            # tract lengths
            samplewise_lengths = all_spectra['homwise_info'][:, 0]
            samplewise_ints = [np.asarray(spect, dtype=np.int32) for spect in samplewise_lengths]
        return samplewise_tract_endpoints, samplewise_ints

    def calc_overlap_size(match, test):
        m1, m2 = match
        t1, t2 = test
        if m1 >= t1 and m2 > t2:  # falls off the rhs
            return t2 - m1
        if m1 < t1 and m2 <= t2:  # falls off the lhs
            return m2 - t1
        if m1 >= t1 and m2 <= t2:  # contained within
            return m2 - m1
        if m1 < t1 and m2 > t2:  # falls off both sides
            return t2 - t1

    def calc_overlap_info(overlap, this_tract):
        length = calc_overlap_size(overlap[:2], this_tract)
        recrate = overlap[2]
        return length, recrate

    def find_tract_recrate(this_tract):
        overlap_infos = [calc_overlap_info(overlap, this_tract)
                         for overlap in recrate_interlap.find(this_tract)]
        lengths, recrates = zip(*overlap_infos)
        avg_recrate = np.average(recrates, weights=lengths)
        return avg_recrate

    def count_tracts_by_bins(one_hom_spectrum, bins):
        bindex = np.digitize(one_hom_spectrum, bins=bins)
        one_hom_binned_counts = np.asarray([len(one_hom_spectrum[bindex == i])
                                            for i in range(num_bins)])
        return one_hom_binned_counts

    this_df = subset_df(df, c=chromosome, d=dominance)
    binned_counts = []
    num_replicates = 0
    all_lengths = []
    all_recrates = []

    for f in this_df["filename"]:
        endpoints, spectrum = import_samplewise_info(f)

        # count the tracts
        these_binned_counts = [count_tracts_by_bins(hom, bins) for hom in spectrum]
        binned_counts.append(these_binned_counts)
        print(f'{num_replicates}: Finished counting tracts.')

        # calculate the recombination rate within each coverage tract
        rep_endpoints = sum(endpoints, [])
        rep_lengths = np.concatenate(spectrum).ravel()
        assert len(rep_endpoints) == len(rep_lengths)
        rep_recrates = [find_tract_recrate(this_tract) for this_tract in rep_endpoints]
        print(f'{num_replicates}: Finished calculating recrates.')

        all_lengths.append(rep_lengths)
        all_recrates.append(rep_recrates)
        num_replicates += 1

    binned_counts = np.reshape(np.asarray(binned_counts), (-1, num_bins))

    return np.sum(binned_counts, axis=0), binned_counts, (all_lengths, all_recrates)


a_counts, _, a_len_rec = count_samplewise_tracts_with_recrates(df_no_sexbias,
                                                               '1', bins,
                                                               recrate_interlap)
x_counts, _, x_len_rec = count_samplewise_tracts_with_recrates(df_no_sexbias,
                                                               'X', bins,
                                                               recrate_interlap)

if save_fig4_data:
    with open('fig4_data/fig4_tracts_A.pkl', 'wb') as file:
        pickle.dump(x_ax, file)
        pickle.dump((a_counts, x_counts), file)
        pickle.dump((exp_cdfA, exp_cdfX), file)
    with open('fig4_data/fig4_tracts_B.pkl', 'wb') as file:
        pickle.dump(a_len_rec, file)
        pickle.dump(x_len_rec, file)


# %% Panel C

# Calculate cumulative proportion of simulated coverage tracts in length bins

x_props = x_counts / sum(x_counts)
x_props = np.cumsum(x_props)
x_props = x_props[1:]

a_props = a_counts / sum(a_counts)
a_props = np.cumsum(a_props)
a_props = a_props[1:]

if save_fig4_data:
    with open('fig4_data/fig4_tracts_C.pkl', 'wb') as file:
        pickle.dump((a_props, x_props), file)


# %% Panel D

def import_spectrum(filepath, kind='introg_spectrum'):
    with np.load(filepath) as all_spectra:
        introg_spectrum = all_spectra[kind]
    return introg_spectrum


def three_spectra_by_replicate(df, dominance, malefrac, autosome='1'):
    dfX = subset_df(df, d=dominance, c='X', m=malefrac)
    dfA = subset_df(df, d=dominance, c=autosome, m=malefrac)
    spectX = [import_spectrum(f) for f in dfX["filename"]]
    spectA = [import_spectrum(f) for f in dfA["filename"]]
    spectrum = np.append(np.concatenate(spectX), np.concatenate(spectA))
    return spectA, spectX, spectrum


def get_cutoffs_by_rep(df, dominance, malefrac, threshold=0.95, norm=1):
    spectsA, spectsX, _ = three_spectra_by_replicate(df, dominance, malefrac)
    cutoffsA = np.asarray([np.quantile(sp, threshold) for sp in spectsA])
    cutoffsX = np.asarray([np.quantile(sp, threshold) for sp in spectsX])
    return cutoffsA / norm, cutoffsX / norm


cmfs = np.asarray([0, 0.5, 1])
boxplot_cutoffs_add = [get_cutoffs_by_rep(df, '0.5SLiM', mf, threshold=0.95)
                       for mf in cmfs]

if save_fig4_data:
    with open('fig4_data/fig4_tracts_D.pkl', 'wb') as file:
        pickle.dump(boxplot_cutoffs_add, file)
