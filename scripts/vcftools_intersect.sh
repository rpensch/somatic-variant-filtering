#!/bin/bash

spm_vcf=$1
sim_vcf=$2
spm_bed=$3
sim_bed=$4
spm_out=$5
sim_out=$6

vcftools --gzvcf $spm_vcf --bed $spm_bed --recode --recode-INFO-all --stdout | gzip -c > $spm_out
vcftools --gzvcf $sim_vcf --bed $sim_bed --recode --recode-INFO-all --stdout | gzip -c > $sim_out
