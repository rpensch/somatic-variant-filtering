import pandas as pd

configfile: "config/config.yaml"

samples = list(pd.read_table(config["samples"]).set_index("sample", drop=False)["sample"])

rule mutect2_filtering:
    input:
        "1_mutect2/raw/Mutect2_filtered_{sample}.vcf.gz"
    output:
        spm="1_mutect2/filtered/{sample}.mutect2.spm.filtered.vcf.gz",
        sim="1_mutect2/filtered/{sample}.mutect2.sim.filtered.vcf.gz"
    log: 
        "logs/mutect2_filtering/{sample}.log"
    shell:
        "scripts/somatic_variant_filtering.py -f {input} -o {output.spm},{output.sim} --spim_separate"

rule strelka_filtering:
    input:
        spm="2_strelka/raw/StrelkaBP_{sample}_somatic_snvs.vcf.gz",
        sim="2_strelka/raw/StrelkaBP_{sample}_somatic_indels.vcf.gz"
    output:
        spm="2_strelka/filtered/{sample}.strelka.spm.filtered.vcf.gz",
        sim="2_strelka/filtered/{sample}.strelka.sim.filtered.vcf.gz"
    log: 
        "logs/strelka_filtering/{sample}.log"
    shell:
        "scripts/somatic_variant_filtering.py -f {input.spm},{input.sim} -o {output.spm},{output.sim} --spim_separate"

rule strelka2bed:
    input:
        spm="2_strelka/filtered/{sample}.strelka.spm.filtered.vcf.gz",
        sim="2_strelka/filtered/{sample}.strelka.sim.filtered.vcf.gz"
    output:
        spm="2_strelka/filtered_bed/{sample}.strelka.spm.filtered.sorted.bed",
        sim="2_strelka/filtered_bed/{sample}.strelka.sim.filtered.sorted.bed"
    log: 
        "logs/strelka2bed/{sample}.log"
    shell:
        "scripts/vcf2bed.sh {input.spm} {input.sim} {output.spm} {output.sim}"

rule intersect:
    input:
        spm_vcf="1_mutect2/filtered/{sample}.mutect2.spm.filtered.vcf.gz",
        sim_vcf="1_mutect2/filtered/{sample}.mutect2.sim.filtered.vcf.gz",
        spm_bed="2_strelka/filtered_bed/{sample}.strelka.spm.filtered.sorted.bed",
        sim_bed="2_strelka/filtered_bed/{sample}.strelka.sim.filtered.sorted.bed" 
    output: 
        spm="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.spm.filtered.intersect.vcf.gz",
        sim="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.sim.filtered.intersect.vcf.gz"
    log:
        "logs/intersect/{sample}.log"
    shell:
        "scripts/vcftools_intersect.sh {input.spm_vcf} {input.sim_vcf} {input.spm_bed} {input.sim_bed} {output.spm} {output.sim}"

rule germline_filtering:
    input:
        spm="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.spm.filtered.intersect.vcf.gz",
        sim="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.sim.filtered.intersect.vcf.gz",
        germl_spm = conf['germline_data']['spm'],
        germl_sim = conf['germline_data']['sim']
    output:
        spm="4_germline_filtered/{sample}.mutect2.strelka.spm.filtered.intersect.germl.vcf.gz",
        sim="4_germline_filtered/{sample}.mutect2.strelka.sim.filtered.intersect.germl.vcf.gz"
    shadow: "full"
    shell:
        """
        vcftools --gzvcf {input.spm} --bed {input.germl_spm} --recode --recode-INFO-all --stdout | gzip -c > {output.spm}
        vcftools --gzvcf {input.sim} --bed {input.germl_sim} --recode --recode-INFO-all --stdout | gzip -c > {output.sim}
        """

rule sample_summary:
    input:
        m2_raw="1_mutect2/raw/Mutect2_filtered_{sample}.vcf.gz",
        m2_filt_spm="1_mutect2/filtered/{sample}.mutect2.spm.filtered.vcf.gz",
        m2_filt_sim="1_mutect2/filtered/{sample}.mutect2.sim.filtered.vcf.gz",
        st_raw_spm="2_strelka/raw/StrelkaBP_{sample}_somatic_snvs.vcf.gz",
        st_raw_sim="2_strelka/raw/StrelkaBP_{sample}_somatic_indels.vcf.gz",
        st_filt_spm="2_strelka/filtered/{sample}.strelka.spm.filtered.vcf.gz",
        st_filt_sim="2_strelka/filtered/{sample}.strelka.sim.filtered.vcf.gz",
        intersect_spm="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.spm.filtered.intersect.vcf.gz",
        intersect_sim="3_mutect2_strelka_intersect/{sample}.mutect2.strelka.sim.filtered.intersect.vcf.gz",
        germl_spm="4_germline_filtered/{sample}.mutect2.strelka.spm.filtered.intersect.germl.vcf.gz",
        germl_sim="4_germline_filtered/{sample}.mutect2.strelka.sim.filtered.intersect.germl.vcf.gz"
    output:
        "summary/samples/{sample}.vcf_summary.tsv"
    shell:
        """
        scripts/sample_summary.vcf_filtering.py --sample {wildcards.sample} \
        --m2_raw {input.m2_raw} --m2_filt {input.m2_filt_spm},{input.m2_filt_sim} \
        --st_raw {input.st_raw_spm},{input.st_raw_sim} --st_filt {input.st_filt_spm},{input.st_filt_sim} \
        --intersect {input.intersect_spm},{input.intersect_sim} \
        --germl {input.germl_spm},{input.germl_sim} -o {output}
        """

rule summary:
    input:
        expand("summary/samples/{sample}.vcf_summary.tsv", sample=samples)
    output:
        "summary/filtering.vcf_summary.tsv"
    shell:
        "scripts/total_summary.vcf_filtering.py {output} {input}"

rule all:
    input:
        "summary/filtering.vcf_summary.tsv"