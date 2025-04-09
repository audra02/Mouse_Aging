import os
import pandas as pd

configfile: "config/config.yaml"
samples = pd.read_table("config/samples.tsv").set_index("sample", drop=False)
GENOME = os.path.abspath("data/reference/bowtie2_index/mm10")

rule all:
    input:
        expand("results/peaks/{sample}_peaks.narrowPeak", sample=samples.index),
        "results/multiqc_report.html"

rule fastqc:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "read1"]
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip = "results/fastqc/{sample}_fastqc.zip"
    conda: "envs/chipseq.yaml"
    threads: 4
    shell:
        "fastqc -t {threads} {input.r1} --outdir results/fastqc"

rule trim_galore:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "read1"]
    output:
        trimmed = "results/trimmed/{sample}_trimmed.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    conda: "envs/chipseq.yaml"
    threads: 4
    shell:
        "trim_galore --gzip -o results/trimmed {input.r1} > {log} 2>&1"

rule bowtie2:
    input:
        r1 = "results/trimmed/{sample}_trimmed.fq.gz"
    output:
        bam = "results/aligned/{sample}.bam"
    params:
        opts = config["params"]["alignment"]["bowtie2"]
    log:
        "logs/bowtie2/{sample}.log"
    shell:
        "bowtie2 {params.opts} -p {threads} -x {GENOME} "
        "-U {input.r1} 2> {log} | "
        "samtools view -bS - | samtools sort -o {output.bam}"

rule macs2:
    input:
        bam = "results/aligned/{sample}.bam"
    output:
        peaks = "results/peaks/{sample}_peaks.narrowPeak",
        summits = "results/peaks/{sample}_summits.bed"
    params:
        opts = config["params"]["chipseq"]["macs2"]
    log:
        "logs/macs2/{sample}.log"
    shell:
        "macs2 callpeak -t {input.bam} -f BAM -g mm -n {wildcards.sample} "
        "--outdir results/peaks {params.opts} > {log} 2>&1"

rule multiqc:
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=samples.index),
        expand("logs/bowtie2/{sample}.log", sample=samples.index),
        expand("logs/macs2/{sample}.log", sample=samples.index)
    output:
        "results/multiqc_report.html"
    conda: "envs/chipseq.yaml"
    shell:
        "multiqc results/ -o results/"
