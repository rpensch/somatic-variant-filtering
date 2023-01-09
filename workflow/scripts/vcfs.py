#!/usr/bin/env python3

import gzip
import pandas as pd

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

def count_spims(vcf):
     
    ''' Count the number of spms and sims '''
     
    n_spm = vcf[(vcf['ALT'].str.len()==1) & 
                (vcf['REF'].str.len()==1)].shape[0]
     
    n_sim = vcf[(vcf['ALT'].str.len()!=1) | 
                (vcf['REF'].str.len()!=1)].shape[0]
     
    return n_spm, n_sim 