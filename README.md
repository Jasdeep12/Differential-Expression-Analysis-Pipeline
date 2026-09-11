# RNA-seq Analysis Pipeline

**This is an educational project meant for my own learning, therefore the outputs and design might not be of fully professional design.

A RNA-seq analysis workflow for paired-end sequence data, built using Python, Snakemake, and Conda.

## Overview

This project takes paired-end RNA-seq reads from *Escherichia coli* K-12 MG1655 through QC, genome alignment, BAM processing, and gene level quantification.

This workflow is designed such that it's reproducible and scalable to multiple samples.

## Workflow

Input : Fastq
-> Fastp
-> HISAT2
-> SAMtools sorting
-> BAM indexing
-> SAMtools flagstat
-> featureCounts
-> MultiQC

## Tools

Snakemake | Workflow management
Fastp | Read Quality Control and Trimming
HISAT2 | Read alignment
SAMtools | BAM sorting and indexing
featureCounts | Gene-level quantification
MultiQC | QC report aggregation
Ncbi-tools-cli | datasets, testing and validation

## Data

Organism:

*Escherichia coli* K-12 MG1655

Example sequencing dataset:

Accession Number: SRR13970441

Reads:
Paired-end FASTQ

## Environment

This pipeline uses a Conda environment defined in:

`envs/RNASeqPipelineProject.yml`

Snakemake directly manages the Conda environments, so no need to activate the environment.

## Example Result

For SRR13970441, HISAT2 achieved a 96.6% overall alignment rate.

The pipeline generated:

- Fastp quality-control reports
- MultiQC summary report
- Sorted and indexed BAM
- SAMtools alignment QC
- featureCounts gene-level count table


## Project Structure

```text
RNASeqPipelineProject/
|- Snakefile
|- README.md
|- config/
|	- samples.tsv
|- data/
|	- raw/
|- envs/
|	- RNASeqPipelineProject.yml
|- reference/
|	- genome.fna
|	- genomic.gff
|	- hisat2_index/
|- results/
|	- bam/
|	- counts/
|	- fastp/
|	- multiqc/
|-scripts/

