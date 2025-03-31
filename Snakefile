# Snakefile

import pandas as pd
import os

# Load config and samples
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"])
samples = samples_df["tissue_age_description"].tolist()

# Auto-detect sample groups
chipseq_samples = [s for s in samples if "K4me3" in s]
rnaseq_samples = [s for s in samples if "RNA" in s]

# Reference files
GENOME_FASTA = config["reference"]["fasta"]
GTF = config["reference"]["gtf"]

# --- Main Workflow ---
rule all:
    input:
        # ChIP-seq outputs
        expand("results/chipseq/peaks/{sample}_peaks.narrowPeak", sample=chipseq_samples),
        expand("results/chipseq/bigwig/{sample}.bw", sample=chipseq_samples),
        
        # RNA-seq outputs
        expand("results/rnaseq/counts/{sample}.counts", sample=rnaseq_samples),
        "results/multiqc/multiqc_report.html"

# --- Shared Rules ---
rule download_sra:
    """
    Downloads SRA files and converts to FASTQ.
    Handles both ChIP-seq and RNA-seq samples.
    """
    output:
        "data/raw/{sample}.fastq.gz"
    params:
        srx = lambda wc: samples_df.loc[samples_df["tissue_age_description"] == wc.sample, "srx_accession"].iloc[0]
    log:
        "logs/download/{sample}.log"
    threads: 2
    conda:
        "envs/sra.yaml"
    shell:
        """
        prefetch {params.srx} -O tmp/ &&
        fasterq-dump tmp/{params.srx} \
            --outdir data/raw/ \
            --skip-technical \
            --threads {threads} \
            -o {wildcards.sample}.fastq &&
        pigz data/raw/{wildcards.sample}.fastq &&
        rm -rf tmp/{params.srx}
        """

rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        html = "results/qc/{sample}_fastqc.html",
        zip = "results/qc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    conda:
        "envs/qc.yaml"
    shell:
        "fastqc {input} -o results/qc/"

# --- ChIP-seq Pipeline ---
rule chipseq_align:
    input:
        fastq = "data/raw/{sample}.fastq.gz",
        index = f"{config['reference']['bowtie2_index']}.1.bt2"
    output:
        "results/aligned/{sample}.bam"
    log:
        "logs/align/{sample}.log"
    threads: 8
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        bowtie2 -x {config[reference][bowtie2_index]} \
            -U {input.fastq} \
            --threads {threads} 2> {log} \
        | samtools sort -o {output}
        samtools index {output}
        """

rule chipseq_call_peaks:
    input:
        bam = "results/aligned/{sample}.bam"
    output:
        "results/chipseq/peaks/{sample}_peaks.narrowPeak"
    log:
        "logs/macs2/{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        macs2 callpeak -t {input.bam} \
            -g mm -n {wildcards.sample} \
            --outdir results/chipseq/peaks/ \
            {config[chipseq][macs2_params]} > {log} 2>&1
        """

# --- RNA-seq Pipeline ---
rule rnaseq_align:
    input:
        fastq = "data/raw/{sample}.fastq.gz",
        index = config["reference"]["star_index"]
    output:
        "results/aligned/{sample}.Aligned.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}.log"
    threads: 8
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN {threads} \
            --outFileNamePrefix results/aligned/{wildcards.sample}. \
            > {log} 2>&1
        """

rule rnaseq_count:
    input:
        bam = "results/aligned/{sample}.Aligned.sortedByCoord.out.bam",
        gtf = GTF
    output:
        "results/rnaseq/counts/{sample}.counts"
    log:
        "logs/counts/{sample}.log"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        featureCounts -a {input.gtf} \
            -o {output} \
            -T 4 \
            {input.bam} > {log} 2>&1
        """

# --- Reporting ---
rule multiqc:
    input:
        expand("results/qc/{sample}_fastqc.zip", sample=samples),
        expand("logs/macs2/{sample}.log", sample=chipseq_samples),
        expand("logs/star/{sample}.log", sample=rnaseq_samples)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc results/ -o results/multiqc/"
