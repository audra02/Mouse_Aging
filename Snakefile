SAMPLES = ["SRX2880744", "SRX2880753", "SRX2880752", "SRX2880736", "SRX2880743", "SRX2880742"]
GENOME_INDEX = "mm10_index"
THREADS = 8

rule all:
    input:
        expand("results/{sample}.counts.txt", sample=SAMPLES)

rule download:
    output:
        "data/{sample}.fastq.gz"
    params:
        sample=lambda wildcards: wildcards.sample
    shell:
        "fastq-dump --split-files --gzip --outdir data {params.sample}"

rule trim:
    input:
        "data/{sample}.fastq.gz"
    output:
        "data/{sample}_trimmed.fastq.gz"
    shell:
        "trim_galore --quality 15 --length 35 --gzip -o data {input}"

rule align:
    input:
        "data/{sample}_trimmed.fastq.gz"
    output:
        "results/{sample}.bam"
    shell:
        "STAR --runThreadN {THREADS} --genomeDir {GENOME_INDEX} --readFilesIn {input} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/{wildcards.sample}"

rule count:
    input:
        "results/{sample}.bam"
    output:
        "results/{sample}.counts.txt"
    shell:
        "featureCounts -T {THREADS} -a annotation.gtf -o results/{wildcards.sample}.counts.txt {input}"

rule analysis:
    input:
        expand("results/{sample}.counts.txt", sample=SAMPLES)
    output:
        "results/differential_expression_results.csv"
    run:
        import pandas as pd
        from scipy.stats import ttest_ind
        
        # Load count data
        dfs = {sample: pd.read_csv(f"results/{sample}.counts.txt", sep='\t', comment='#') for sample in SAMPLES}
        
        # Example differential expression analysis (basic t-test between samples)
        genes = dfs[SAMPLES[0]].iloc[:, 0]
        results = {"gene": genes}
        for i in range(1, len(SAMPLES), 2):
            group1 = dfs[SAMPLES[i-1]].iloc[:, -1]
            group2 = dfs[SAMPLES[i]].iloc[:, -1]
            results[f"{SAMPLES[i-1]}_vs_{SAMPLES[i]}"] = [ttest_ind(group1, group2, equal_var=False).pvalue for _ in genes]
        
        df_results = pd.DataFrame(results)
        df_results.to_csv("results/differential_expression_results.csv", index=False)
