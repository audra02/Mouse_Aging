# Mouse_Aging
## ChIP-seq Analysis Pipeline

Snakemake pipeline for processing H3K4me3 ChIP-seq data (mouse) from raw FASTQs to peak calling.

### Workflow Steps

1. **Quality Control**  
   - `FastQC`: Initial read quality reports  
   - `TrimGalore`: Adapter trimming & quality filtering  
     - Outputs: `{sample}_trimmed.fq.gz`, `{sample}.fastq.gz_trimming_report.txt`

2. **Alignment**  
   - `Bowtie2`: Maps reads to mm10 genome  
   - `samtools`: BAM sorting/indexing  
     - Outputs: `{sample}.bam`, `{sample}.bam.bai`

3. **Peak Calling**  
   - `MACS2`: Identifies H3K4me3 peaks (narrowPeak format)  
     - Outputs: `{sample}_peaks.narrowPeak`, `{sample}_summits.bed`

4. **Reporting**  
   - `MultiQC`: Aggregates all QC metrics

### Configuration
samples:
  - SRR1979921
  - SRR1979925
  - SRR1979927
  - SRR1979930

genome: mm10

Conda environments (see envs/chipseq.yaml)
