#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:08:23 2023

@author: egibson
"""

import pandas as pd

# %% make outgroup .json

sample_info_file = '/Users/egibson/Documents/science/Grad/SP20/demog20/proj/neand-purging/hmmix/1kG_tutorial/integrated_call_samples_v3.20130502.ALL.panel'

sdf = pd.read_table(sample_info_file)

afrdf = sdf[sdf['super_pop'] == 'AFR']

ogdf = sdf[sdf['pop'].isin(['YRI', 'ESN', 'MSL'])]  # LWK and GWD included in rev1

femaledf = ogdf[ogdf['gender'] == 'female']

jsondf = femaledf.set_index('sample')

with open('/Users/egibson/Documents/science/Grad/SP20/demog20/proj/neand-purging/hmmix/1kG_local_attempt/outgroup_rev2_females.json', 'w') as outgroup_file:
    jsondf.to_json(outgroup_file, orient="split", indent=4)

# %% and ingroup .json

female_ingroup_samples = midf[midf['is_female']]
jsondf = female_ingroup_samples.set_index('name')

with open('/Users/egibson/Documents/science/Grad/SP20/demog20/proj/neand-purging/hmmix/1kG_local_attempt/ingroup_rev2_females.json', 'w') as outgroup_file:
    jsondf.to_json(outgroup_file, orient="split", indent=4)
