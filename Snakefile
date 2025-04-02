import os
configfile: "config/config.yaml"

# Sample handling
samples = list(pd.read_csv(config["samples"])["tissue_age_description"])
chipseq_samples = [s for s in samples if "K4me3" in s]
rnaseq_samples = [s for s in samples if "RNA" in s]

rule all:
    input:
        # ChIP-seq outputs
        expand("results/chipseq/peaks/{sample}_peaks.narrowPeak", sample=chipseq_samples),
        # RNA-seq outputs
        expand("results/rnaseq/counts/{sample}.counts", sample=rnaseq_samples),
        # QC
        "results/multiqc_report.html"

# --- Core Analysis Rules ---
rule align_chipseq:
    input:
        fastq = "data/raw/{sample}.fastq.gz",
        index = config["reference"]["bowtie2_index"] + ".1.bt2"
    output:
        bam = "results/aligned/{sample}.bam",
        bai = "results/aligned/{sample}.bam.bai"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120  # minutes
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        bowtie2 -x {config[reference][bowtie2_index]} \
            -U {input.fastq} \
            --threads {threads} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule call_peaks:
    input:
        "results/aligned/{sample}.bam"
    output:
        "results/chipseq/peaks/{sample}_peaks.narrowPeak"
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input} -g mm -n {wildcards.sample} --outdir results/chipseq/peaks/"

rule align_rnaseq:
    input:
        fastq = "data/raw/{sample}.fastq.gz",
        index = config["reference"]["star_index"]
    output:
        "results/aligned/{sample}.Aligned.sortedByCoord.out.bam"
    threads: 12
    resources:
        mem_mb=32000,
        runtime=240
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix results/aligned/{wildcards.sample}.
        """

# --- Utility Rules ---
rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "results/qc/{sample}_fastqc.html"
    resources:
        runtime=30
    conda:
        "envs/qc.yaml"
    shell:
        "fastqc {input} -o results/qc/"

rule multiqc:
    input:
        expand("results/qc/{sample}_fastqc.html", sample=samples)
    output:
        "results/multiqc_report.html"
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc results/ -o results/"
