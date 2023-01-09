#!/usr/bin/env python3

import os.path
import pandas as pd
import argparse
from vcfs import get_vcf_header, load_vcf, count_spims
'''
m2_raw="1_mutect2/raw/Mutect2_filtered_T-BlueSkye_vs_N-BlueSkye.vcf.gz"
m2_filt="1_mutect2/filtered/T-BlueSkye_vs_N-BlueSkye.mutect2.spm.filtered.vcf.gz,1_mutect2/filtered/T-BlueSkye_vs_N-BlueSkye.mutect2.sim.filtered.vcf.gz"
st_raw="2_strelka/raw/StrelkaBP_T-BlueSkye_vs_N-BlueSkye_somatic_snvs.vcf.gz,2_strelka/raw/StrelkaBP_T-BlueSkye_vs_N-BlueSkye_somatic_indels.vcf.gz"
st_filt="2_strelka/filtered/T-BlueSkye_vs_N-BlueSkye.strelka.spm.filtered.vcf.gz,2_strelka/filtered/T-BlueSkye_vs_N-BlueSkye.strelka.sim.filtered.vcf.gz"
intersect="3_mutect2_strelka_intersect/T-BlueSkye_vs_N-BlueSkye.mutect2.strelka.spm.filtered.intersect.vcf.gz,3_mutect2_strelka_intersect/T-BlueSkye_vs_N-BlueSkye.mutect2.strelka.sim.filtered.intersect.vcf.gz"
steps = [(m2_raw,'Mutect2_raw'), 
         (m2_filt,'Mutect2_filtered'), 
         (st_raw,'Strelka_raw'), 
         (st_filt, 'Strelka_filtered'), 
         (intersect, 'Mutect2_Strelka_intersect')]


'''

# Parse command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True, type=str)
parser.add_argument('--m2_raw', required=True, type=str)
parser.add_argument('--m2_filt', required=True, type=str)
parser.add_argument('--st_raw', required=True, type=str)
parser.add_argument('--st_filt', required=True, type=str)
parser.add_argument('--intersect', required=True, type=str)
parser.add_argument('--germl', required=True, type=str)
parser.add_argument('-o', '--out', required=True, type=str)

args = parser.parse_args()

steps = [(args.m2_raw,'Mutect2_raw'), 
         (args.m2_filt,'Mutect2_filtered'), 
         (args.st_raw,'Strelka_raw'), 
         (args.st_filt, 'Strelka_filtered'), 
         (args.intersect, 'Mutect2_Strelka_intersect'),
         (args.germl, 'germline_filtered')]

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