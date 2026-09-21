# Differential Expression Analysis Pipeline

A Differential Expression analysis workflow for paired-end and single-end sequence data, built using Python, R, Snakemake, and Conda.

## Overview

This project takes single-ended and paired-end RNA-seq reads from a sample of your choice through QC, genome alignment, BAM processing, and gene level quantification, and DESeq2 analysis.

This workflow is designed such that it's reproducible and scalable to multiple samples and conditions.

## Commands

Samples are registered in `config/samples.tsv` using `add_sample.py`, a Click-based CLI with 4 subcommands.

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

### condition

Adds a separate column to the samples.tsv file for control(s) and treatment(s)/condition(s), this column allows the DESeq2 rule to automatically match and perform analysis.

```bash
python add_sample.py condition --control 3 --treatment 3 --treatment 3 --names penicillin --names ampicillin
```
The above command will add control to the condition column for the first 3 rows, and penicillin to rows 4-6, and ampicillin to rows 7-9.

```bash
python add_sample.py condition --control 2 --treatment 2
```
The above command will add control to the condition column for the first 2 rows, and will add treatment1 to rows 3-4.


## Workflow

Input : Fastq
-> Fastp
-> HISAT2
-> SAMtools sorting
-> BAM indexing
-> SAMtools flagstat
-> featureCounts
-> MultiQC
-> DESeq2

## Usage
Can be run with the following command:
```bash
snakemake --cores 4 --use-conda
```
The HISAT2 index is not provided nor is the annotation file.

## Tools

Snakemake | Workflow management
Fastp | Read Quality Control and Trimming
HISAT2 | Read alignment
SAMtools | BAM sorting and indexing
featureCounts | Gene-level quantification
MultiQC | QC report aggregation
DESeq2 | Differential Expression Analysis
Ncbi-tools-cli | datasets, testing and validation

## Data

Originally tested on (will no longer be sufficient due to DESeq2):

Organism:

*Escherichia coli* K-12 MG1655

Original sequencing dataset:

Accession Number: SRR13970441

Recommended Testing Set:

https://www.ncbi.nlm.nih.gov/bioproject/PRJEB75208
Runs 9-14

Reads:
Single-end FASTQ,
Paired-end FASTQ

## Environment

This pipeline uses Conda environments defined in:

`envs/RNASeqPipelineProject.yml`
`envs/deseq2.yml`


Snakemake directly manages the Conda environments, so no need to activate the environment as long as you have Snakemake installed.


## Example Result

For SRR13970441, HISAT2 achieved a 96.6% overall alignment rate.

The pipeline generated:

- Fastp quality-control reports
- MultiQC summary report
- Sorted and indexed BAM
- SAMtools alignment QC
- featureCounts gene-level count table
- Count Matrix for input into DESeq2
- DESeq2 Analysis


## Project Structure

```text
RNASeqPipelineProject/
|- Snakefile
|- add_sample.py
|- README.md
|- config/
|	- config.yaml
|	- samples.tsv
|- data/
|	- raw/
|- envs/
|	- RNASeqPipelineProject.yml
|	- deseq2.yml
|- reference/
|	- hisat2_index/
|- results/
|	- bam/
|	- counts/
|	- fastp/
|	- multiqc/
|-scripts/
|	- DE.R
|	- make_matrix.py
```
## Citations

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Chen, J. W., Shrestha, L., Green, G., Leier, A., & Marquez-Lago, T. T. (2023). The hitchhikers' guide to RNA sequencing and functional analysis. Briefings in bioinformatics, 24(1), bbac529. https://doi.org/10.1093/bib/bbac529

Genovese, M., De Santis, M., Giordano, A., Marotta, P., & Galietta, L. J. V. (2026). TRPV4-driven ATP release activates CFTR through ADORA2B and P2RY2 receptors in the airway epithelium. iScience, 29(9), 117304. https://doi.org/10.1016/j.isci.2026.117304

