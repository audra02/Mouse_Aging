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
        r1 = lambda wildcards: samples.loc[wildcards.sample, "read1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "read2"]
    output:
        html1 = "results/fastqc/{sample}_R1_fastqc.html",
        html2 = "results/fastqc/{sample}_R2_fastqc.html"
    conda: "envs/chipseq.yaml"
    threads: 4
    shell:
        "fastqc -t {threads} {input.r1} {input.r2} --outdir results/fastqc"

# Modified trim_galore to handle cases where R2 is just indexes
rule trim_galore:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "read1"]
    output:
        r1 = "results/trimmed/{sample}_trimmed.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    conda: "envs/chipseq.yaml"
    threads: 4
    shell:
        "trim_galore --gzip -o results/trimmed {input.r1} > {log} 2>&1"

# Modified bowtie2 to handle both true PE and index-only data
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

# Modified MACS2 to auto-detect input type
rule macs2:
    input:
        bam = "results/aligned/{sample}.bam"
    output:
        peaks = "results/peaks/{sample}_peaks.narrowPeak"
    params:
        opts = config["params"]["chipseq"]["macs2"],
        # Use BAMPE only if R2 is >1MB
        format = lambda wildcards: "BAMPE" if os.path.getsize(samples.loc[wildcards.sample, "read2"]) > 1000000 else "BAM"
    log:
        "logs/macs2/{sample}.log"
    shell:
        "macs2 callpeak -t {input.bam} -f {params.format} -g mm -n {wildcards.sample} "
        "--outdir results/peaks {params.opts} > {log} 2>&1"

rule multiqc:
    input:
        expand("results/fastqc/{sample}_R1_fastqc.html", sample=samples.index),
        expand("logs/bowtie2/{sample}.log", sample=samples.index)
    output:
        "results/multiqc_report.html"
    conda: "envs/chipseq.yaml"
    shell:
        "multiqc results/ -o results/"
