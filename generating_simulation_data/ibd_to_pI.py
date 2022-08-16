#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 20:58:16 2020

@author: egibson
"""

import tskit as tsk
import numpy as np
import sys
import pickle

# sys.path.append('/path/to/ms-repo/generating_simulation_data/')
import tree_tools as tt


# modified from https://stackoverflow.com/a/15273749
def union_tracts(a):
    b = []
    for begin, end in sorted(a):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    return b


def calculate_pI_using_IBD(slim_ts_file, recomb_file, hom_popn=2,
                           min_tract_length=None, max_extant_samples=None,
                           save_ibd_dicts=False, sex=None,
                           timecourse_calc=False, throw_muts=False):
    # Load tree sequence
    slim_ts = tsk.load(slim_ts_file)
    seq_length = slim_ts.sequence_length

    # Determine time of introgression, oldest Neanderthals (for max IBD time),
    #   and ancestral pop'n size (for recapitation) from samples.
    sample_times = [slim_ts.node(s).time for s in slim_ts.samples()]
    del slim_ts
    sample_times, sample_counts = np.unique(sample_times, return_counts=True)

    if timecourse_calc:
        mixtime = sample_times[-3]
    else:
        mixtime = sample_times[1]
    oldest_neands = mixtime + 1
    Nanc = sample_counts[-1]  # in haplotypes, not indivs
    if sex == 'X':
        Nanc = round(Nanc * 0.75)

    print("Successfully read from original ts.  Moving to pyslim.")
    # Recapitate treeseq from SLiM
    ts = tt.process_treeseq(slim_ts_file, recomb_file, initial_Ne=Nanc,
                            throw_muts=throw_muts, dump_recap=False, sex=sex)
    print("Finished recapitation.  Finding simplification samples.")

    # Simplify ts to introgressing Neanderthals and extant samples.
    #   Write out as '.trees.sts'
    if timecourse_calc:
        sts, neanderthals, extant_samples, timecourse_times, timecourse_sample_list = tt.find_neands_and_sample_timecourse(ts, hom_popn, mixtime)
    else:
        sts, neanderthals, extant_samples = tt.find_neanderthal_introg_class(ts,
                                                                             hom_popn,
                                                                             mixtime)
    print(f"Found {len(neanderthals)} introgressing Neanderthals.")
    sts.dump(slim_ts_file.removesuffix('.decap') + '.sts')

    # Restrict number of extant samples in which to calculate pI
    if (max_extant_samples is not None) and (len(extant_samples) > max_extant_samples):
        hom_samples = np.random.choice(extant_samples, size=max_extant_samples,
                                       replace=False)
    else:
        hom_samples = extant_samples

    # Helper function for single sample
    def calculate_pI_within_one_hom(this_hom, return_ibd_dict):
        # Calculate IBD between this_hom and all Neanderthals
        neand_hom_pairs = [(this_hom, n) for n in neanderthals]
        neand_ibd_this_hom = sts.tables.find_ibd(neand_hom_pairs,
                                                 max_time=oldest_neands,
                                                 min_length=min_tract_length)
        endpoints_per_neand = [[(l, r) for l, r in zip(v['left'], v['right'])]
                               for k, v in neand_ibd_this_hom.items()]
        endpoints = sum(endpoints_per_neand, [])  # flatten list of endpoint pairs

        # Find introgression tract spectrum
        introg_tracts = union_tracts(endpoints)
        introg_tract_lengths = [r - l for [l, r] in introg_tracts]
        introg_this_hom = sum(introg_tract_lengths)  # in bp

        pI_this_hom = introg_this_hom / seq_length

        if return_ibd_dict:
            return pI_this_hom, neand_ibd_this_hom
        else:
            return pI_this_hom

    # Calculate pI for each hom sample.
    if save_ibd_dicts:  # Write each sample IBD dictionary into an '.ibd' file
        ibd_megadict = {}
        pI_hom_samples = []
        print("Calculating pI and saving IBD dictionaries.")
        num_homs = len(hom_samples)
        for i, h in enumerate(hom_samples):
            pI, pdict = calculate_pI_within_one_hom(h, return_ibd_dict=True)
            ibd_megadict.update(pdict)
            pI_hom_samples.append(pI)
            if (i < 5) or (i % 25 == 0):
                print(f"\t{'{:4d}'.format(i+1)}/{num_homs} extant samples finished.")
        print(f"\t{'{:4d}'.format(num_homs)}/{num_homs} extant samples finished.")

        outfn = slim_ts_file.removesuffix('.trees.decap') + '.ibd'
        with open(outfn, 'wb') as outf:
            pickle.dump(ibd_megadict, outf)
        print(f"Find IBD file at {outfn}.")

    else:  # Return only the pI estimate for each sample.
        print("Calculating pI.")
        pI_hom_samples = [calculate_pI_within_one_hom(h, return_ibd_dict=False)
                          for h in hom_samples]

    # Write pI for each hom sample to text file
    outfn = slim_ts_file.removesuffix('.trees.decap') + '.pI'
    with open(outfn, 'a') as outf:
        for hom, pI in zip(, pI_hom_samples):
            line = str(int(hom)) + '\t' + str(float(pI))
            outf.write(line)
            outf.write('\n')

    if timecourse_calc:
        # Calculate pI for timecourse samples and write to file
        # Doesn't save IBD dicts.
        # Each line is [time hom pI]
        outfn = slim_ts_file.removesuffix('.trees.decap') + '.timecourse.pI'
        with open(outfn, 'a') as outf:
            gen = 1  # iterator for printing progress
            for time, sample_list in zip(timecourse_times, timecourse_sample_list):
                for hom in sample_list:
                    pI = calculate_pI_within_one_hom(hom, return_ibd_dict=False)

                    line = str(int(time)) + '\t' + str(int(hom)) + '\t' + str(float(pI))
                    outf.write(line)
                    outf.write('\n')
                print(f"\t Timecourse pI calculated for {gen} generations after introgression.")
                gen += 1  # !!! Assumes generation-by-generation remembering.

    return pI_hom_samples


# %% Set parameters, calculate pI, save pI to file, [optional save IBD to file]

def main():
    # Set parameters, taking '.trees' file from sys
    ts_file = str(sys.argv[1])
    recomb_file = str(sys.argv[2])
    sex = str(sys.argv[3])  # 'X' for chrX, 'A' for aut
    timecourse_calc = str(sys.argv[4])


    min_tract_length = 500  # (bp) per Steinrucken
    max_extant_samples = 1000  # Restrict number of samples in which to calc pI.
    save_ibd_dicts = True
    throw_muts = False

    hom_popn = 2

    # Calculate pI for each sample and write to '.pI' text file
    calculate_pI_using_IBD(ts_file, recomb_file, hom_popn=hom_popn,
                           min_tract_length=min_tract_length,
                           max_extant_samples=max_extant_samples,
                           save_ibd_dicts=save_ibd_dicts, sex=sex,
                           timecourse_calc=timecourse_calc,
                           throw_muts=throw_muts)


if __name__ == "__main__":
    main()
