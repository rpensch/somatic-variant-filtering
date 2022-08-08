#!/bin/bash

spm=$1
sim=$2
spm_out=$3
sim_out=$4

# Convert spm vcfs to bed files
printf "#CHROM\tSTART\tSTOP\n" > $spm_out
zcat $spm | convert2bed -i vcf --snvs | 
cut -f1-3 | sort -k1,1 -k2,2n \
>> $spm_out

# Convert sim vcfs to bed files
printf "#CHROM\tSTART\tSTOP\n" > $sim_out
zcat $sim | convert2bed -i vcf --insertions | 
cut -f1-3 | sort -k1,1 -k2,2n \
>> $sim_out
zcat $sim | convert2bed -i vcf --deletions | 
cut -f1-3 | sort -k1,1 -k2,2n \
>> $sim_out
