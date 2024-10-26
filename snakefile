# Snakefile

import os
import glob

# Set the working directory
working_dir = "/home/akshay/Akshay/RNAseq/Snakemake/"
os.chdir(working_dir)

# Define input and output files
FASTQ_FILES = glob.glob("*.fastq.gz")
TRIMMED_R1 = [f"trimmed_{r1}" for r1 in glob.glob("*_R1.fastq")]
TRIMMED_R2 = [f"trimmed_{r2}" for r2 in glob.glob("*_R2.fastq")]
SAM_FILES = [f"{r1.replace('_R1.fastq', '.sam')}" for r1 in TRIMMED_R1]
BAM_FILES = [f"{sam_file.replace('.sam', '_sorted.bam')}" for sam_file in SAM_FILES]

# Rule to generate MD5 checksums
rule md5_checksums:
    output:
        "md5_checksums.txt"
    shell:
        "md5sum {input} > {output}".format(input=' '.join(FASTQ_FILES))

# Rule to unzip FASTQ files
rule unzip_fastq:
    input:
        FASTQ_FILES
    output:
        expand("{sample}.fastq", sample=[f.replace(".gz", "") for f in FASTQ_FILES])
    shell:
        "gunzip {input}"

# Rule to run FastQC
rule fastqc:
    input:
        glob_wildcards("*.fastq").fastq
    output:
        directory("QC_reports")
    shell:
        "fastqc -o QC_reports/ {input}"

# Rule to run MultiQC
rule multiqc:
    input:
        directory("QC_reports")
    output:
        "QC_reports/multiqc_report.html"
    shell:
        "multiqc -o QC_reports/ QC_reports/"

# Rule to trim paired-end reads using Fastp
rule trim_reads:
    input:
        r1="{sample}_R1.fastq",
        r2="{sample}_R2.fastq"
    output:
        output_r1="trimmed_{sample}_R1.fastq",
        output_r2="trimmed_{sample}_R2.fastq",
        html_report="QC_reports/fastp_{sample}.html",
        json_report="QC_reports/fastp_{sample}.json"
    threads:
        8
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.output_r1} -O {output.output_r2} "
        "-h {output.html_report} -j {output.json_report} --thread {threads}"

# Rule to build Hisat2 index
rule hisat2_index:
    input:
        "reference_genome.fa"
    output:
        "index_prefix.1.ht2"
    shell:
        "hisat2-build {input} index_prefix"

# Rule to align reads with Hisat2
rule hisat2_align:
    input:
        r1="trimmed_{sample}_R1.fastq",
        r2="trimmed_{sample}_R2.fastq"
    output:
        "{sample}.sam"
    shell:
        "hisat2 -x index_prefix -1 {input.r1} -2 {input.r2} -S {output}"

# Rule to convert SAM to BAM, sort, and generate statistics
rule sam_to_bam:
    input:
        "{sample}.sam"
    output:
        "{sample}_sorted.bam"
    shell:
        "samtools view -bS {input} | samtools sort -o {output} && "
        "samtools flagstat {output} > {wildcards.sample}_alignment_statistics.txt"

# Rule to run featureCounts
rule feature_counts:
    input:
        expand("{sample}_sorted.bam", sample=glob_wildcards("*.bam"))
    output:
        "gene_counts.txt"
    shell:
        "featureCounts -T 8 -a annotation.gtf -o {output} -g gene_id -t exon -s 1 -p {input}"

# Specify the order of execution
rule all:
    input:
        "md5_checksums.txt",
        "QC_reports/multiqc_report.html",
        "gene_counts.txt"
