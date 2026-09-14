# RNA-seq Analysis Pipeline

**This is an educational project meant for my own learning, therefore the outputs and design might not be of fully professional design.

A RNA-seq analysis workflow for paired-end sequence data, built using Python, Snakemake, and Conda.

## Overview

This project takes paired-end RNA-seq reads from *Escherichia coli* K-12 MG1655 through QC, genome alignment, BAM processing, and gene level quantification.

This workflow is designed such that it's reproducible and scalable to multiple samples.

## Commands

Samples are registered in `config/samples.tsv` using `add_sample.py`, a Click-based CLI with three subcommands.

### fetch

Downloads a single accession via `prefetch` + `fasterq-dump` and registers it as a sample.

```bash
python add_sample.py fetch SRR13970441
python add_sample.py fetch SRR13970441 --name ecoli_ctrl_1
python add_sample.py fetch SRR13970441 --force
```

- `--name` — sample name to use in `samples.tsv` (defaults to the accession)
- `--force` — overwrite an existing sample with the same name

### fetchall

Downloads two or more accessions in one call. Failures (a failed download, missing paired output, or a duplicate sample name) are skipped rather than stopping the batch, and a summary is printed at the end.

```bash
python add_sample.py fetchall SRR000001 SRR000002 SRR000003
python add_sample.py fetchall SRR000001 SRR000002 --name ctrl_1 --name treat_1
```

- `--name` — one name per accession, given in the same order (defaults to each accession)
- `--force` — overwrite existing samples with matching names

For a single accession, use `fetch` instead.

### local

Registers FASTQ files that are already downloaded, either as individual files or a directory. R1/R2 mates are paired automatically by filename (`_1`/`_2` or `_R1`/`_R2`).

```bash
python add_sample.py local data/raw/mysample_R1.fastq.gz data/raw/mysample_R2.fastq.gz
python add_sample.py local data/raw/
```

- `--force` — overwrite existing samples with matching names


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

