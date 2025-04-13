SAMPLES = ["SRR1979921", "SRR1979925", "SRR1979927", "SRR1979930"]
GENOME = "/scratch/lustre/home/mekl0332/Mouse_Aging/data/reference/bowtie2_index/mm10"

rule all:
    input:
        # FastQC
        expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        # Trimmed files
        expand("results/trimmed/{sample}_trimmed.fq.gz", sample=SAMPLES),
        expand("results/trimmed/{sample}.fastq.gz_trimming_report.txt", sample=SAMPLES),
        # Alignment
        expand("results/aligned/{sample}.bam", sample=SAMPLES),
        expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),
        # Peaks
        expand("results/peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),
        expand("results/peaks/{sample}_summits.bed", sample=SAMPLES),
        # MultiQC
        "results/multiqc_report.html"

# FastQC
rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"  # Removed _R1
    output:
        html="results/fastqc/{sample}_fastqc.html",
        zip="results/fastqc/{sample}_fastqc.zip"
    threads: 4
    shell:
        "mkdir -p results/fastqc && fastqc -t {threads} {input} --outdir results/fastqc"

# TrimGalore
rule trim_galore:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        trimmed=temp("results/trimmed/{sample}_trimmed.fq.gz"),
        report="results/trimmed/{sample}.fastq.gz_trimming_report.txt"
    log:
        "logs/snakemake/trim_galore/{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p results/trimmed logs/snakemake/trim_galore && \
        trim_galore --gzip --basename {wildcards.sample} \
          --output_dir results/trimmed {input} > {log} 2>&1
        """

# Bowtie2 alignment
rule bowtie2:
    input:
        trimmed=rules.trim_galore.output.trimmed
    output:
        bam="results/aligned/{sample}.bam",
        bai="results/aligned/{sample}.bam.bai"
    params:
        opts="--very-sensitive --no-mixed --no-discordant"
    log:
        "logs/snakemake/bowtie2/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p results/aligned logs/snakemake/bowtie2 && \
        bowtie2 {params.opts} -p {threads} -x {GENOME} \
          -U {input.trimmed} 2>> {log} | \
        samtools view -@ 2 -bS - | \
        samtools sort -@ 2 -o {output.bam} - && \
        samtools index {output.bam}
        """

# MACS2
rule macs2:
    input:
        bam=rules.bowtie2.output.bam
    output:
        peaks="results/peaks/{sample}_peaks.narrowPeak",
        summits="results/peaks/{sample}_summits.bed"
    params:
        opts="--qvalue 0.05 --call-summits"
    log:
        "logs/snakemake/macs2/{sample}.log"
    shell:
        """
        mkdir -p results/peaks logs/snakemake/macs2 && \
        macs2 callpeak -t {input.bam} -f BAM -g mm -n {wildcards.sample} \
          --outdir results/peaks {params.opts} > {log} 2>&1
        """

# MultiQC
rule multiqc:
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=SAMPLES),  # Removed _R1
        expand("results/trimmed/{sample}.fastq.gz_trimming_report.txt", sample=SAMPLES),
        expand("logs/snakemake/bowtie2/{sample}.log", sample=SAMPLES),
        expand("logs/snakemake/macs2/{sample}.log", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    shell:
        "multiqc results/ -o results/"
