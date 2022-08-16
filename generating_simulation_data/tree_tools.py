#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import msprime
import numpy as np


def process_treeseq(tree_file, region_info_file, throw_muts=False,
                    sex=None, unif_recomb=False, n_scale=1, initial_Ne=10000,
                    verbose=False, dump_recap=False, mut_rate=1.5e-8,
                    neu_or_neg=0):
    """
    Recapitates .trees output from SLiM simulation and further processes:
        Unless simulating the Xchr, we overlay neutral mutations using msprime.
        If simulating the Xchr, we simplify to remove all Y chromosomes.
    Overwrites .trees file UNLESS tree_file has extension '.orig'

    Parameters
    ----------
    tree_file : str
        Path to .trees file output from SLiM simulation.
    region_info_file : str
        Path to "sim_seq_info_[region].txt" file with recRate info.
    neu_or_neg : int, optional
        Null model of neutral variation or alt. model of rec. del. var.
        neu = 2; neg = 0;  per "dominance" parameter in original code.
        The default is 0.
    n_scale : int, optional
        Scaling factor. The default is 10.
    initial_Ne : int, optional
        Initial population size of founding population (p1).
        The default is 10000.
    verbose : bool, optional
        Print info about treeseqs before and after recapitation.
        The default is False.

    Output
    ------
    Writes out .trees file, recapitated, with EITHER neutral mutations
    overlayed OR Y chr removed.  Currently OVERWRITES the original .trees file
    UNLESS the original has ext '.orig'.

    Returns
    ------
    ts : treeSeq

    """
    # Load ts
    slim_ts = pyslim.load(tree_file, legacy_metadata=True)
    print("ts loaded into pyslim.")

    # Set recombination map
    def make_region_recombination_map(region_filename):
        if unif_recomb:
            recomb_map = msprime.RateMap(unif_recomb)
        else:
            positions = []
            rates = []
            with open(region_filename, 'r') as file:
                for line in file:
                    if "recRate" in line:
                        components = line.split(" ")
                        positions.append(int(components[1]))
                        rates.append(float(components[2]))
            # adapted from https://pyslim.readthedocs.io/en/latest/tutorial.html#recapitation-with-a-nonuniform-recombination-map
            # step 1
            positions.insert(0, 0)
            # step 2
            # rates.append(0.0)
            # step 3
            positions[-1] += 1
            recomb_map = msprime.RateMap(position=positions, rate=rates)
        return recomb_map

    # recomb_map = make_region_recombination_map(region_info_file)

    # Recapitate
    # original ts loaded with pyslim, so a SLiMTreeSeq object.
    # https://pyslim.readthedocs.io/en/latest/python_api.html#pyslim.SlimTreeSequence.recapitate
    n_p1 = round(initial_Ne / n_scale)
    # print("Beginning recapitation with recombination map.")
    recap_ts = slim_ts.recapitate(recombination_map=recomb_map, Ne=n_p1)
    assert max([t.num_roots for t in recap_ts.trees()]) == 1

    if sex == 'X':
        # simplify to remove Y chromosomes
        chrXs = [recap_ts.node(s).id for s in recap_ts.samples()
                 if recap_ts.node(s).metadata['genome_type'] == pyslim.GENOME_TYPE_X]
        ts = recap_ts.simplify(samples=chrXs)
    else:
        ts = recap_ts

    if throw_muts:
        # Manually create msprime mutation times by editing tables ( :o )
        tables = ts.dump_tables()
        final_slim_generation = ts.slim_generation
        time_agos = [final_slim_generation - one_tmut.metadata['mutation_list'][0]['slim_time']
                     for one_tmut in tables.mutations]
        tables.mutations.time = np.array(time_agos)
        ts = tables.tree_sequence()

        # Need mutation rate variation to ensure equal density of mutations along chromosome
        # mutation rate scaling notes from SLiM sim code:
        #   // exon uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
        # 	// to achieve an overall mutation rate of 1.5e-8, need 2.31/(1+2.31) fraction of all muts to be nonsyn
        # 	// i.e. ~0.6979 fraction of all muts should be the deleterious ones simulated here as "m1".
        # 	// the remaining ~0.3021 fraction of all mut.s should be thrown later as neutrals.
        # HOWEVER, this logic applies if we were only simulating exon regions.
        def make_region_mutation_map(region_filename, exon_rate, intron_rate):
            starts = []
            ends = []
            with open(region_filename, 'r') as file:
                for line in file:
                    if "exon" in line:
                        components = line.split(" ")
                        starts.append(int(components[1]))
                        ends.append(int(components[2]))
            ends = [e + 1 for e in ends]  # starts of introns
            positions = [val for pair in zip(starts, ends) for val in pair]
            # https://stackoverflow.com/a/7946825
            positions.insert(0, 0)
            positions.append(slim_ts.sequence_length)
            rates = [intron_rate, exon_rate] * len(starts)
            rates.append(intron_rate)

            mut_map = msprime.RateMap(position=positions, rate=rates)
            return mut_map

        exon_mut_rate = (1 / (1 + 2.31)) * mut_rate * n_scale  # already sim'd rest as deleterious
        intron_mut_rate = mut_rate * n_scale  # no muts in intronic regions yet

        mut_rate_map = make_region_mutation_map(region_info_file, exon_mut_rate,
                                                intron_mut_rate)

        def throw_mutations(ts, rate_obj, from_slim_gen=None, to_slim_gen=None,
                            mut_type=5):
            # Creating a SLiM mutation model for msprime
            max_id = -1
            for mut in ts.mutations():
                for d in mut.derived_state.split(","):
                    max_id = max(max_id, int(d))
            model = msprime.SLiMMutationModel(type=mut_type, next_id=max_id + 1)

            # working out timings
            if to_slim_gen is not None:
                start_time = slim_ts.slim_generation - to_slim_gen
            else:
                start_time = None
            if from_slim_gen is not None:
                end_time = slim_ts.slim_generation - from_slim_gen
            else:
                end_time = None

            # Throw mutations with msprime 1
            ts = msprime.sim_mutations(ts, model=model, rate=rate_obj,
                                       start_time=start_time, end_time=end_time,
                                       keep=True)
            return ts

        ts = throw_mutations(ts, mut_rate, to_slim_gen=0)
        ts = throw_mutations(ts, mut_rate_map, from_slim_gen=0, mut_type=6)

    if dump_recap:
        # Write out .trees file
        # default is to overwrite...
        out_name = tree_file
        # ... unless SLiM ts file has extension .decap
        if tree_file[-6:] == '.decap':
            out_name = tree_file[:-6]
        ts.dump(out_name)

    # Messing around
    if verbose:
        import matplotlib.pyplot as plt
        # Mess around with original ts
        print(f"There are {slim_ts.num_individuals} individuals in the SLiM tree.")
        founders = slim_ts.first_generation_individuals()
        initial_pops = [slim_ts.individual(f).population for f in founders]
        print(f"Population(s) represented in the founders: {np.unique(initial_pops)}.")
        # numpy array of how long ago each individual was born
        # SLiM automatically remembers the individuals that comprise the first generation of any new subpopulation created with addSubpop(), for easy recapitation and other analysis.
        indiv_times = slim_ts.individual_times
        plt.figure()
        plt.hist(indiv_times)
        recorded_times = np.unique(indiv_times)
        print(f"When indivs. were added to the ts, in time ago: {recorded_times}.")
        root_num_spectrum = [t.num_roots for t in slim_ts.trees()]
        # plt.figure()
        # plt.hist(root_num_spectrum)
        max_roots = max(root_num_spectrum)
        print(f"Before recapitation, the max number of roots was {max_roots}.")
        # Mess around with full recap'd, mutated ts
        recap_root_num_spectrum = [t.num_roots for t in ts.trees()]
        recap_root = max(recap_root_num_spectrum)
        print(f"After recapitation, the max number of roots is {recap_root}.")
        print("As for throwing mutations...")
        print(f"before there were {slim_ts.num_mutations}; after, {ts.num_mutations}. ")

    pyts = pyslim.SlimTreeSequence(ts, legacy_metadata=True)
    return pyts


def translate_from_slim_pop(tree_seq, s_ids=None):
    popkey = [(p.id, p.metadata.slim_id) for p in tree_seq.populations()
              if p.metadata is not None]
    pops = [p for (p, s) in popkey]
    slims = [s for (p, s) in popkey]
    if s_ids is None:  # just get all populations with a slim_id
        s_ids = slims
    indices_of_slim_pops = np.where(np.isin(slims, s_ids))[0]
    return tuple(np.asarray(pops)[indices_of_slim_pops])


def simplify_since_admixture(ts, recip_popn, time_since_adm):
    # check that we got the recapitated treeseq
    assert max(t.num_roots for t in ts.trees()) == 1

    recip_popn = translate_from_slim_pop(ts, s_ids=[recip_popn])

    samples_to_keep = [ts.node(s).id for s in ts.samples(population=recip_popn)
                       if ts.node(s).time <= time_since_adm]

    sts = ts.simplify(samples=samples_to_keep)
    print("Finished simplification.")
    return sts


def find_neanderthal_introg_class(ts, recip_popn, time_since_adm, precomputed_sts=None):
    if precomputed_sts is not None:
        sts = pyslim.load(precomputed_sts)
    else:
        sts = simplify_since_admixture(ts, recip_popn, time_since_adm)

    neanderthals = [sts.node(s).id for s in sts.samples() if sts.node(s).time == time_since_adm]
    extant_samples = [sts.node(s).id for s in sts.samples() if sts.node(s).time == 0]

    return sts, neanderthals, extant_samples  # these are neanderthals and extant homs


def find_neands_and_sample_timecourse(ts, recip_popn, time_since_adm,
                                      max_timecourse_sample_number=1000):
    sts = simplify_since_admixture(ts, recip_popn, time_since_adm)

    neanderthals = [sts.node(s).id for s in sts.samples() if sts.node(s).time == time_since_adm]
    extant_samples = [sts.node(s).id for s in sts.samples() if sts.node(s).time == 0]

    samples = sts.samples()
    times = np.asarray([sts.node(s).time for s in sts.samples()])
    sample_times = np.unique(times)

    timecourse_times = sample_times[1:-1]  # remove introg. and exants
    list_of_timecourse_samples = []
    for gen in timecourse_times:
        these_samples = samples[np.where(times == gen)]
        if len(these_samples) > max_timecourse_sample_number:
            these_samples = np.random.choice(these_samples,
                                             size=max_timecourse_sample_number,
                                             replace=False)
        list_of_timecourse_samples.append(these_samples)

    return sts, neanderthals, extant_samples, timecourse_times, list_of_timecourse_samples
