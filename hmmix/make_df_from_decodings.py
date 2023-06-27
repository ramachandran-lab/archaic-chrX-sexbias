#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:47:19 2023

@author: egibson
"""
import pandas as pd
from interlap import InterLap

# Filter tracts
restrict_to_archaic_SNPs = True  # this applies to autosomes too!  all the others are just for chrX atm.  also does not change denominator
remove_PAR = True
remove_selection = True
remove_centromere = False
remove_gaps = False

output_file_basename = 'output.csv.gz'


# %%  Define functions
def process_chrX(df, remove_PAR=False, remove_selection=False, remove_centromere=False, remove_gaps=False):

    def filter_regions(df, filter_starts, filter_ends):
        inter = InterLap()
        par_endpoints = list(zip(filter_starts, filter_ends))
        inter.update(par_endpoints)
        # determine which archaic tracts overlap
        df['endpoints'] = list(zip(df['start'], df['end']))
        not_under_mask = [coords not in inter for coords in df['endpoints']]
        # frac_chrX_tracts_under_mask = 1 - sum(not_under_mask) / len(not_under_mask)
        tot_bp_chrX_before_filtering = df['length'].sum()
        # drop all tracts that touch a filter region
        df = df.loc[not_under_mask]
        tot_bp_chrX_after_filtering = df['length'].sum()
        tot_bp_under_mask = tot_bp_chrX_before_filtering - tot_bp_chrX_after_filtering
        pct_filtered_this_step = tot_bp_under_mask / tot_bp_chrX_before_filtering * 100
        print(f"\tremoved {pct_filtered_this_step}% of chrX coverage (bp)")
        return df, pct_filtered_this_step

    if remove_centromere:
        centromere_start = [58632012]
        centromere_end = [61632012]
        print("filtering centromere:")
        df, _ = filter_regions(df, centromere_start, centromere_end)

    if remove_gaps:
        # all gaps from https://genome.ucsc.edu/cgi-bin/hgTables.  group: All Tables, table: gap
        gaps_file = 'all_hg19_GRCh37_chrX_gaps.csv'
        gdf = pd.read_csv(gaps_file, sep='\t')
        gap_starts = gdf['chromStart']
        gap_ends = gdf['chromEnd']
        print("filtering gaps:")
        df, _ = filter_regions(df, gap_starts, gap_ends)

    if remove_PAR:
        # starts and ends from https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/par.txt
        PAR_starts = [60001, 154931044]
        PAR_ends = [2699520, 155260560]
        print("filtering PAR:")
        df, _ = filter_regions(df, PAR_starts, PAR_ends)

    if remove_selection:
        selection_file = 'skov_selection_targets.csv'
        sdf = pd.read_csv(selection_file, skiprows=1, skip_blank_lines=True)
        sdf = sdf.dropna()
        selection_starts = sdf['Start.2'] * 100000
        selection_ends = sdf['End.2'] * 100000
        print("filtering selection")
        df, _ = filter_regions(df, selection_starts, selection_ends)

    return df


def find_segments_on_one_hap(indiv_ID, hap, chr_type='autosome',
                             remove_PAR=False, remove_selection=False, remove_centromere=False, remove_gaps=False,
                             restrict_to_archaic_SNPs=False):

    decode_output_file = f'decoded_{chr_type}s/decoded.{indiv_ID}.{chr_type}.{hap}.txt'

    df = pd.read_table(decode_output_file)

    #  retain only archaic segments with >= 80% support
    df = df.loc[(df['state'] == 'Archaic') & (df['mean_prob'] >= 0.8)]

    #  remove segments with no archaic variants
    if restrict_to_archaic_SNPs:
        df = df[df['admixpopvariants'] > 0]

    # further process PAR and selection targets
    if (chr_type == 'chrX'):
        df = process_chrX(df, remove_PAR=remove_PAR, remove_selection=remove_selection, remove_centromere=remove_centromere, remove_gaps=remove_gaps)

    df['hap'] = hap
    df['name'] = indiv_ID
    return df


def find_segments_in_one_indiv(indiv_ID, chr_type='autosome', remove_PAR=False,
                               remove_selection=False, remove_centromere=False, remove_gaps=False, restrict_to_archaic_SNPs=False):
    return pd.concat([find_segments_on_one_hap(indiv_ID, hap, chr_type=chr_type,
                                               remove_PAR=remove_PAR,
                                               remove_selection=remove_selection, remove_centromere=remove_centromere, remove_gaps=remove_gaps,
                                               restrict_to_archaic_SNPs=restrict_to_archaic_SNPs)
                      for hap in ('hap1', 'hap2')])


def summarize_length_in_one_indiv(indiv_ID, chr_type='autosome', remove_PAR=False,
                                  remove_selection=False, remove_centromere=False, remove_gaps=False, restrict_to_archaic_SNPs=False):
    df = find_segments_in_one_indiv(indiv_ID, chr_type=chr_type,
                                    remove_PAR=remove_PAR,
                                    remove_selection=remove_selection, remove_centromere=remove_centromere, remove_gaps=remove_gaps,
                                    restrict_to_archaic_SNPs=restrict_to_archaic_SNPs)
    return df.groupby(['name']).sum(numeric_only=True)

# indiv_df.groupby(['name', 'hap']).sum()


def summarize_length(all_indiv_IDs, chr_type='autosome', remove_PAR=False,
                     remove_selection=False, remove_centromere=False, remove_gaps=False, restrict_to_archaic_SNPs=False):
    return pd.concat([summarize_length_in_one_indiv(indiv_ID, chr_type=chr_type,
                                                    remove_PAR=remove_PAR,
                                                    remove_selection=remove_selection, remove_centromere=remove_centromere, remove_gaps=remove_gaps,
                                                    restrict_to_archaic_SNPs=restrict_to_archaic_SNPs)
                      for indiv_ID in all_indiv_IDs])


# %%  summarize tract length per individual
id_file = 'ingroup_haps.txt'
id_list = [line.rstrip() for line in open(id_file)]

aut_length_df = summarize_length(id_list)
chrX_length_df = summarize_length(id_list, chr_type='chrX',
                                  remove_PAR=remove_PAR,
                                  remove_selection=remove_selection,
                                  remove_centromere=remove_centromere, 
                                  remove_gaps=remove_gaps,
                                  restrict_to_archaic_SNPs=restrict_to_archaic_SNPs)


# %%  merge aut and chrX info

aut_length_df.rename(columns={"length": "aut_len"}, inplace=True)
chrX_length_df.rename(columns={"length": "x_len"}, inplace=True)
df = pd.concat([aut_length_df, chrX_length_df], axis=1)

# %% add sample info, including sex

sample_info_file = 'integrated_call_samples_v3.20130502.ALL.panel'

sdf = pd.read_table(sample_info_file, index_col=0, usecols=(0, 1, 2, 3),
                    names=('name', 'pop', 'superpop', 'sex'), skiprows=1)
sdf['is_female'] = sdf['sex'].apply(lambda x: x == 'female')

df = df.join(sdf)


# %% calculate per base pair rates and ratio
# Get lengths from Hg19 (aka GRCh37) to make per-bp numbers
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
aut_len_hg19 = 2981033286
x_len_hg19 = 155270560

PAR_length = 2969035
selection_length = 13300000
centromere_length = 3000000
gaps_length = 4170000  # gdf['size'].sum()

output_file_prefix = ''
x_len_to_remove = 0

if restrict_to_archaic_SNPs:
    output_file_prefix += 'yesarchaicvar_'
if remove_PAR:
    x_len_to_remove += PAR_length
    output_file_prefix += 'noPAR_'
if remove_selection:
    x_len_to_remove += selection_length
    output_file_prefix += 'noselection_'
if remove_centromere:
    x_len_to_remove += centromere_length
    output_file_prefix += 'nocentromere_'
if remove_gaps:
    x_len_to_remove += gaps_length
    output_file_prefix += 'nogaps_'

output_name = output_file_prefix + output_file_basename
x_len_tot = x_len_hg19 - x_len_to_remove


df['aut_pbp'] = df['aut_len'] / (2 * aut_len_hg19)
df['x_pbp'] = df['x_len'] / (2 * x_len_tot)
df['aut_x_ratio'] = df['aut_pbp'] / df['x_pbp']


# %% write out as .csv

df.to_csv(output_name, index_label='name', na_rep='NaN')

