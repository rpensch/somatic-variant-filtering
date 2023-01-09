#!/usr/bin/env python3

import os.path
import pandas as pd
import argparse
from vcfs import get_vcf_header, load_vcf, count_spims

# Parse command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True, type=str)
parser.add_argument('--m2_raw', required=True, type=str)
parser.add_argument('--m2_filt', required=True, type=str)
parser.add_argument('--st_raw', required=True, type=str)
parser.add_argument('--st_filt', required=True, type=str)
parser.add_argument('--intersect', required=True, type=str)
parser.add_argument('--germl', required=False, type=str)
parser.add_argument('-o', '--out', required=True, type=str)

args = parser.parse_args()

if args.germl:
    germl_args = [(args.germl, 'germline_filtered')]
else:
    germl_args = []

steps = [(args.m2_raw,'Mutect2_raw'), 
         (args.m2_filt,'Mutect2_filtered'), 
         (args.st_raw,'Strelka_raw'), 
         (args.st_filt, 'Strelka_filtered'), 
         (args.intersect, 'Mutect2_Strelka_intersect')] + germl_args

spim_columns = [f for flatten in list(zip(
               [step[1] + '_spm' for step in steps], 
               [step[1] + '_sim' for step in steps])) for f in flatten]

summary = pd.DataFrame(columns=spim_columns, index=[args.sample])

for s in steps:
    
    # Load vcf file(s)
    vcf = pd.DataFrame(columns = get_vcf_header(s[0].split(',')[0]))
    for path in s[0].split(','):
        if os.path.exists(path):
            vcf = vcf.append(load_vcf(path))
        else:
            raise IOError(f'{path} cannot be found. Did you provide a correct path?')
    
    n_spm, n_sim = count_spims(vcf)
    summary.loc[args.sample][s[1]+'_spm'] = n_spm
    summary.loc[args.sample][s[1]+'_sim'] = n_sim

summary.to_csv(args.out, sep='\t', index_label='sample')