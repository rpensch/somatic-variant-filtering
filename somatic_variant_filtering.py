#!/usr/bin/env python3

import os.path
import pandas as pd
import argparse
import gzip

def get_vcf_header(vcf_path):
     
    ''' Get the vcf header '''
     
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:  
            if line.startswith("#CHROM"):
                vcf_names = [n.strip('\n') for n in line.split('\t')]
                break
    return vcf_names

def save_vcf_comments(vcf_path, out_path):
     
    ''' Get the vcf comment rows '''
     
    with gzip.open(vcf_path, 'rt') as f:
        comments = []
        for line in f:
            if line.startswith("##"):
                comments.append(line)       
            if line.lower().startswith("#CHROM"):
                break
                    
    with gzip.open(out_path, 'wt') as o:
        for line in comments:
            o.write(line)

def load_vcf(vcf_path):
         
    ''' Load vcf file '''
     
    vcf_names = get_vcf_header(vcf_path)
    vcf = pd.read_csv(vcf_path, comment='#', sep='\t', header=None, names=vcf_names)
     
    return vcf

def filter_vcf(vcf):
     
    ''' Remove low quality calls and duplicates '''
     
    filtered = vcf[(vcf['FILTER'].str.lower()=='pass') | 
                   (vcf['FILTER'].str.lower()=='.')] \
                    .drop_duplicates(['#CHROM','POS','REF','ALT'])
     
    return filtered

def check_sim_ratio(vcf, max):
     
    ''' Check if the sim ratio exeeds max '''
     
    n_spm = vcf[(vcf['ALT'].str.len()==1) & 
                (vcf['REF'].str.len()==1)].shape[0]
     
    n_sim = vcf[(vcf['ALT'].str.len()!=1) | 
                (vcf['REF'].str.len()!=1)].shape[0]
     
    ratio = n_sim / (n_spm + n_sim)
     
    if ratio > max:
        import warnings
        warnings.warn(f"WARNING: Somatic indel ratio exceeds {max}.")

def separate_spims(vcf):
     
    ''' Split variants into somatic point and indel mutations '''
     
    spm = vcf[(vcf['ALT'].str.len()==1) & 
              (vcf['REF'].str.len()==1)]
     
    sim = vcf[(vcf['ALT'].str.len()!=1) | 
              (vcf['REF'].str.len()!=1)]
     
    return spm, sim    

# Parse command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', required=True, type=str)
parser.add_argument('-o', '--out', required=True, type=str)
parser.add_argument('--max_sim', default=0.2, type=float)
parser.add_argument('--spim_separate', action='store_true')

args = parser.parse_args()

file = args.file
out = args.out

if args.spim_separate:
    assert len(args.out.split(',')) == 2 , "If using --spim_separate, provide two output files separated by ',' (spm first)."

# Load vcf file(s)
vcf = pd.DataFrame(columns = get_vcf_header(file.split(',')[0]))
for path in file.split(','):
    if os.path.exists(path):
        vcf = vcf.append(load_vcf(path))
    else:
        raise IOError(f'{path} cannot be found. Did you provide a correct path?')

# Filter vcfs and check sim ratio
filtered = filter_vcf(vcf)
check_sim_ratio(filtered, args.max_sim)

# If spim_separate, save spms and sims separately, else all
if args.spim_separate:

    spm, sim = separate_spims(filtered)

    for i, spim in enumerate([spm, sim]):

        if len(file.split(',')) == 2:
            j = i
        elif len(file.split(',')) == 1:
            j = 0

        save_vcf_comments(file.split(',')[j], out.split(',')[i])
        spim.to_csv(out.split(',')[i], sep='\t', index=False, mode='a')

else:
    save_vcf_comments(file.split(',')[0], out.split(',')[0])
    filtered.to_csv(out.split(',')[0], sep='\t', index=False, mode='a')

