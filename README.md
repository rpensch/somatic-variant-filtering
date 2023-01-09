# Somatic variant filtering pipeline

This repository contains a simple snakemake pipeline to filter somatic point mutations (SPMs) and somatic indel mutations (SIMs). It uses raw variant vcf files from two variant calles (Mutect2 and Strelka), removes low quality as well as duplicate variants and intersects the results from both callers to obtain only high-confidence variant calls. Optionally, a germline resource can be used to remove known common germline variants. 

Pipeline steps:

1. Filtering

In a first step, variants that do not pass the quality filters applied by the callers are removed. Occasionally, callers will list a variant twice for the same sample which can be a problem for some downstream analyses. These duplicates are removed as well. Finally, Mutect2 variants are separated into SPMs and SIMs (Strelka vcfs are separated by default). A warning will be raised if the percentage of indels of the total variants of a sample exeeds 20 %. 

2. Intersection

Mutect2 and Strelka filtered variants are intersected and only variants that are found by both callers are kept. 

3. Optional: Germline filter

Common germline variants from a given bed file can be filtered out in this step.

4. Summary

Summary statistics with variant counts for each sample at each step are calculated.

## Installation

For instructions on how to set up snakemake, see [snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).

Clone this repository:

    `git clone https://github.com/rpensch/somatic_variant_filtering.git`

## Required tools

The following tools are required to run somatic variant filtering with this pipelie:

- snakemake (tested version 7.18.2)
- python3 (tested version 3.9.5)
- vcftools (tested version 0.1.16)
- BEDOPS (tested version 2.4.39)

## Running the pipeline

### 1. Input

To run the pipeline provide a tsv file specifying sample names and vcf file paths with the following columns (including header line):

- sample
- mutect2_vcf
- strelka_snv_vcf
- strelka_indel_vcf

Optional: bed files of common germline variants to be excluded

### 2. Config

The config file can be found at `config/config.yaml`.

In the config file:

1. Specify the location of the input file in the config file under `input_files`. 
2. If you want to run the pipeline with the optional germline filter:
    - Set `germline_filter` to `true` and
    - Provide bed files for germline snps and indels under `germline_data`, `snps` and `indels`

### 3. Run

Make sure you are in the somatic_variant_filtering directory:

    cd somatic_variant_filtering/

Test if the set up was successful with a dry run:

    snakemake -np all

Run the pipeline with 1 core:

    snakemake --cores 1 all

### 4. Output

Output files can be found in the results directory. The final filtered vcf files are in `results/3_mutect2_strelka_intersect` and separate for somatic point mutations(SPMs) and somatic indel mutations (SIMs).

- mysample.mutect2.strelka.spm.filtered.intersect.vcf.gz
- mysample.mutect2.strelka.sim.filtered.intersect.vcf.gz

If a germline filter was used, the results will be found in `results/4_germline_filtered/`.

The `results/summary` directory contains a summary of variant counts for all samples.