#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 13:31:29 2023

@author: egibson
"""
import sys
import pandas as pd

from hmm_functions import DecodeModel, read_HMM_parameters_from_file
from helper_functions import Load_observations_weights_mutrates

from collections import defaultdict
import numpy as np
from numba import njit
import json
import math

from helper_functions import find_runs, Annotate_with_ref_genome, Make_folder_if_not_exists, flatten_list
# -----------------------------------------------------------------------------
# Running on an individual from 1000 genomes
# -----------------------------------------------------------------------------


def convert_params_to_chrX(indiv_ID, ploidy='haploid'):
    if ploidy == 'haploid':
        parameter_file = f'trained_files/trained.{indiv_ID}.haploid.json'
        params = read_HMM_parameters_from_file(parameter_file)
        params.emissions *= (5 / 6)
    elif ploidy == 'diploid':
        parameter_file = f'diploid_trained_files/trained.{indiv_ID}.diploid.json'
        params = read_HMM_parameters_from_file(parameter_file)

        # determine sex of individuals
        sample_info_file = 'integrated_call_samples_v3.20130502.ALL.panel'
        sdf = pd.read_table(sample_info_file, index_col=0, usecols=(0, 1, 2, 3),
                            names=('name', 'pop', 'superpop', 'sex'), skiprows=1)
        is_female = sdf['sex'].apply(lambda x: x == 'female')

        if is_female[indiv_ID]:
            params.emissions *= (5 / 6)
        else:
            params.emissions[1] = (5 / 6) * (params.emissions[1] - (params.emissions[0] / 2))
            params.emissions[0] = (5 / 6) * (params.emissions[0] / 2)
    return params


def Write_Decoded_output(outputprefix, segments, obs_file = None, admixpop_file = None, extrainfo = False):

    # # Load archaic data
    # if admixpop_file is not None:
    #     admix_pop_variants, admixpop_names = Annotate_with_ref_genome(admixpop_file, obs_file)

    # Are we doing haploid/diploid?
    outfile_mapper = {}
    for _, _, _, _, _, _, _, ploidity, _ in segments:
        if outputprefix == '/dev/stdout':
            outfile_mapper[ploidity] = '/dev/stdout'
        else:
            outfile_mapper[ploidity] = f'{outputprefix}.{ploidity}.txt'

    # Make output files and write headers
    outputfiles_handlers = defaultdict(str)
    for ploidity, output in outfile_mapper.items():

        Make_folder_if_not_exists(output)
        outputfiles_handlers[ploidity] = open(output, 'w')
        out = outputfiles_handlers[ploidity]

        # if admixpop_file is not None:
        #     if extrainfo:
        #         out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tadmixpopvariants\t{}\tvariant\tfoundin\n'.format('\t'.join(admixpop_names)))
        #     else:
        #         out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tadmixpopvariants\t{}\n'.format('\t'.join(admixpop_names)))
        # else:
        out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\n')

    # Go through segments and write to output
    for chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, ploidity, variants in segments:
        out = outputfiles_handlers[ploidity]

        # if admixpop_file is not None:
        #     archiac_variants_dict = defaultdict(int)
        #     for snp_position in variants.split(','):
        #         variant = admix_pop_variants[f'{chrom}_{snp_position}']
        #         if variant != '':
        #             if '|' in variant:
        #                 for ind in variant.split('|'):
        #                     archiac_variants_dict[ind] += 1
        #             else:
        #                 archiac_variants_dict[variant] += 1

        #             archiac_variants_dict['total'] += 1

        #     archaic_variants = '\t'.join([str(archiac_variants_dict[x]) for x in ['total'] + admixpop_names])

        #     if extrainfo:
        #         for variant in variants.split(','):

        #             if admix_pop_variants[f'{chrom}_{variant}'] == '':
        #                 foundin = 'none'
        #             else:
        #                 foundin = admix_pop_variants[f'{chrom}_{variant}']

        #             print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, archaic_variants, variant, foundin, sep = '\t', file = out)
        #     else:
        #         print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, archaic_variants, sep = '\t', file = out)

        # else:
        if chrom == 'X':  # added by ETC 20230122
            print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, sep='\t', file=out)

    # Close output files
    for ploidity, out in outputfiles_handlers.items():
        out.close()
    pass


def write_outfile_from_segments(indiv_ID, segments, ploidy='haploid'):
    if ploidy == 'haploid':
        outfile_prefix = f'decoded_chrXs/decoded.{indiv_ID}.chrX'
    else:
        outfile_prefix = f'diploid_decoded_chrXs/decoded.{indiv_ID}.chrX'
    Write_Decoded_output(outfile_prefix, segments)
    pass


def main():
    args = sys.argv[1:]
    indiv_ID = args[0]
    ploidy = args[1]

    is_haploid = (ploidy == 'haploid')

    obs, chroms, starts, variants, mutrates, weights = Load_observations_weights_mutrates(obs_file = f'obs_files/obs.{indiv_ID}.txt',
                                                                                          weights_file = 'strickmask.bed',
                                                                                          mutrates_file = 'mutationrate.bed',
                                                                                          window_size = 1000,
                                                                                          haploid = is_haploid)

    hmm_parameters = convert_params_to_chrX(indiv_ID, ploidy=ploidy)
    segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)
    write_outfile_from_segments(indiv_ID, segments, ploidy=ploidy)
    pass


if __name__ == "__main__":
    main()
