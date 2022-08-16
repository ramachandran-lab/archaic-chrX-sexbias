#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:04:03 2020

@author: egibson
"""

import pandas as pd
import os
import glob
import numpy as np


def get_properties_from_pI_filename(filename):
    # TODO: switch to .removesuffix/prefix for Python 3.9+

    parts = filename.split('_')
    assert len(parts) >= 5

    try:
        assert 'chr' in parts[0]
        chromosome = parts[0].strip('chr')
    except AssertionError:
        print(f"Can't read chromosome from {filename}.")
        chromosome = ''

    try:
        assert 'h' in parts[1]
        dominance = parts[1].strip('h')
    except AssertionError:
        print(f"Can't read dominance from {filename}.")
        dominance = ''

    try:
        assert 'sscale' in parts[2]
        s_scale = int(parts[2].strip('sscale'))
    except AssertionError:
        print(f"Can't read s_scale from {filename}.")
        s_scale = 9999

    try:
        assert 'NsNr' in parts[3]
        NsNr = round(float(parts[3].strip('NsNr')), ndigits=2)
    except AssertionError:
        print(f"Can't read NsNr from {filename}.")
        NsNr = 9999

    # seed is the last part
    seed = int(parts[-1].strip('.timecourse.pI'))

    return chromosome, dominance, s_scale, NsNr, seed


def write_df_for_one_pI_file(filename, path_to_pI_files, include_timecourse=False):

    path_to_pI = path_to_pI_files + filename

    if '.timecourse' in filename:
        if include_timecourse is True:
            df = pd.read_csv(path_to_pI, sep='\t', names=['timepoint', 'node',
                                                          'pI'])
        else:  # ignore this .timepoint.pI file
            return
    else:
        df = pd.read_csv(path_to_pI, sep='\t', names=['node', 'pI'])

    chrom, dom, sscale, nsnr, seed = get_properties_from_pI_filename(filename)

    df.insert(0, "chromosome", chrom)
    df.insert(1, "dominance", dom)
    df.insert(2, "s_scale", sscale)
    df.insert(3, "NsNr", nsnr)
    df.insert(4, "seed", seed)
    df["filename"] = filename

    return df


def compile_df_from_pI_files(path_to_pI_files, include_timecourse=False):
    pI_paths = glob.iglob(path_to_pI_files + '*.pI')
    pI_files = [os.path.basename(p) for p in pI_paths]

    dfs = [write_df_for_one_pI_file(f, path_to_pI_files,
                                    include_timecourse=include_timecourse)
           for f in pI_files]
    pIdf = pd.concat(dfs)

    if include_timecourse is True:
        # Replace missing timepoints with 0
        pIdf['timepoint'] = pIdf['timepoint'].fillna(0)

    return pIdf


def subset_pI(df, c=None, d=None, s=None, n=None, r=None, t=None):
    masks = []
    for val, col in zip((c, d, s, n, r, t),
                        ("chromosome", "dominance", "s_scale", "NsNr", "seed",
                         "timepoint")):
        if val is not None:
            if isinstance(val, list):
                masks.append(df[col].isin(val))
            else:
                masks.append(df[col] == val)

    if len(masks) == 0:  # nothing to subset
        return df
    elif len(masks) == 1:  # just one series to subset on
        return df[masks[0]]
    else:  # subset multiple columns by combining series
        masks = np.vstack(masks)
        megamask = np.all(masks, axis=0)
        return df[megamask]
