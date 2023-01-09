#!/usr/bin/env python3

import os.path
import pandas as pd
import argparse
from vcfs import get_vcf_header, load_vcf, filter_vcf, check_sim_ratio, separate_spims, save_vcf_comments

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