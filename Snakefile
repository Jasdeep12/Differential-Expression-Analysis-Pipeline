import pandas as pd
import re


configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = samples["sample"].tolist()


def is_paired(sample):
	r2 = samples.loc[samples["sample"] == sample, "R2"].iloc[0]
	return pd.notna(r2) and str(r2).strip() != ""

def trimmed_reads(sample):
	if is_paired(sample):
		return [f"results/trimmed/{sample}_1.fastq.gz", f"results/trimmed/{sample}_2.fastq.gz"]
	else:
		return [f"results/trimmed/{sample}.fastq.gz"]
		
def fastp_input(wildcards):
	row = samples.loc[samples["sample"] == wildcards.sample].iloc[0]
	if is_paired(wildcards.sample):
		return {"r1": row["R1"], "r2": row["R2"]}
	return {"r1": row["R1"]}
	
def fastp_output(wildcards):
	if is_paired(wildcards.sample):
		return {
			"trim1": f"results/trimmed/{wildcards.sample}_1.fastq.gz",
			"trim2": f"results/trimmed/{wildcards.sample}_2.fastq.gz",
			"json": f"results/fastp/{wildcards.sample}.json",
			"html": f"results/fastp/{wildcards.sample}.html",
		}
	return {
		"trim1": f"results/trimmed/{wildcards.sample}.fastq.gz",
		"json": f"results/fastp/{wildcards.sample}.json",
		"html": f"results/fastp/{wildcards.sample}.html",
	}


HISAT2_INDEX = config["reference"]["hisat2_index"]
ALIGN_THREADS = config["threads"]["align"]
SAMTOOLS_THREADS = config["threads"]["samtools"]

rule all:
	input:
		[f for sample in SAMPLES for f in trimmed_reads(sample)],
		expand("results/fastp/{sample}.json",
		sample=SAMPLES
	),	
		expand("results/fastp/{sample}.html",
		sample=SAMPLES
	),		
		expand("results/bam/{sample}.sorted.bam.bai",	
		sample=SAMPLES
	),	
		expand("results/qc/{sample}.flagstat.txt",
		sample=SAMPLES
	),
		expand("results/counts/{sample}_counts.txt",
		sample=SAMPLES
	),
		expand("results/counts/{sample}_counts.txt.summary",
		sample=SAMPLES
	),
		"results/multiqc/multiqc_report.html",
		"results/deseq/deseq_results.csv",
		"results/volcano/volcano.pdf",
		"results/volcano/volcano.png"
									
									
PAIRED_SAMPLES = [s for s in SAMPLES if is_paired(s)]
SINGLE_SAMPLES = [s for s in SAMPLES if not is_paired(s)]
paired_regex = "|".join(re.escape(s) for s in PAIRED_SAMPLES) or "(?!)"
single_regex = "|".join(re.escape(s) for s in SINGLE_SAMPLES) or "(?!)"

rule fastp_pe:
	input:
		r1=lambda w: samples.loc[samples["sample"] == w.sample, "R1"].iloc[0],
		r2=lambda w: samples.loc[samples["sample"] == w.sample, "R2"].iloc[0]
	output:
		trim1="results/trimmed/{sample}_1.fastq.gz",
		trim2="results/trimmed/{sample}_2.fastq.gz",
		json="results/fastp/{sample}.json",
		html="results/fastp/{sample}.html"
	wildcard_constraints:
		sample=paired_regex
	conda:
		"envs/RNASeqPipelineProject.yml"
	threads: ALIGN_THREADS
	shell:
		"fastp -i {input.r1} -I {input.r2} -o {output.trim1} -O {output.trim2} "
		"-h {output.html} -j {output.json} --thread {threads}"

rule fastp_se:
	input:
		r1=lambda w: samples.loc[samples["sample"] == w.sample, "R1"].iloc[0]
	output:
		trim1="results/trimmed/{sample}.fastq.gz",
		json="results/fastp/{sample}.json",
		html="results/fastp/{sample}.html"
	wildcard_constraints:
		sample=single_regex
	conda:
		"envs/RNASeqPipelineProject.yml"
	threads: ALIGN_THREADS
	shell:
		"fastp -i {input.r1} -o {output.trim1} "
		"-h {output.html} -j {output.json} --thread {threads}"

def align_input(wildcards):
    reads = trimmed_reads(wildcards.sample) 
    inputs = {
        "index": f"{HISAT2_INDEX}.1.ht2",
        "r1": reads[0]
    }
    
    # Add the second read pair if the sample is paired-end
    if is_paired(wildcards.sample):
        inputs["r2"] = reads[1]
        
    return inputs

rule align:
	input:
		unpack(align_input)
	output:
		bam="results/bam/{sample}.sorted.bam"
	conda:
		"envs/RNASeqPipelineProject.yml"	
	log:
		"logs/{sample}.hisat2.log"
	
	threads: ALIGN_THREADS

	run:
		if is_paired(wildcards.sample):
			shell(
				"hisat2 -p {threads} -x {HISAT2_INDEX} -1 {input.r1} -2 {input.r2} "
				"2> {log} | samtools sort -@ {SAMTOOLS_THREADS} -o {output.bam} -"
			)
		else:
			shell(
				"hisat2 -p {threads} -x {HISAT2_INDEX} -U {input.r1} "
				"2> {log} | samtools sort -@ {SAMTOOLS_THREADS} -o {output.bam} -"
			)


rule index_bam:
	input:
		"results/bam/{sample}.sorted.bam"

	output:
		"results/bam/{sample}.sorted.bam.bai"
	conda:
		"envs/RNASeqPipelineProject.yml"

	shell:
		"""
		samtools index {input} {output}
		"""


rule bam_qc:
	input:
		bam="results/bam/{sample}.sorted.bam"
	
	output:
		"results/qc/{sample}.flagstat.txt"
	
	conda:
		"envs/RNASeqPipelineProject.yml"

	shell:
		"""
		samtools flagstat {input.bam} > {output}
		"""


rule quantify:
	input: 
		bam="results/bam/{sample}.sorted.bam",
		annotation=config["reference"]["annotation"]
	output:
		counts="results/counts/{sample}_counts.txt",
		summary="results/counts/{sample}_counts.txt.summary"
	
	conda:
		"envs/RNASeqPipelineProject.yml"
	params:
		paired_flag = lambda wildcards: "-p" if is_paired(wildcards.sample) else ""
	shell:
		"""
		featureCounts {params.paired_flag} \
			-a {input.annotation} \
			-o {output.counts} \
			-t gene \
			-g gene_id \
			{input.bam}
		"""

rule multiqc:
	input:
		fastpjson=expand("results/fastp/{sample}.json", sample=SAMPLES),
		flagstat=expand("results/qc/{sample}.flagstat.txt", sample=SAMPLES),
		hisat2=expand("logs/{sample}.hisat2.log", sample=SAMPLES)
	output:
		html="results/multiqc/multiqc_report.html"

	conda:
		"envs/RNASeqPipelineProject.yml"

	shell:
		"""
		multiqc \
			results/fastp \
			results/qc \
			logs \
			--outdir results/multiqc \
			--force
		"""

rule matrix:
	input:
		counts=expand("results/counts/{sample}_counts.txt", sample=SAMPLES)
	output:
		matrix="results/matrix/count_matrix.tsv"
	params:
		samples=SAMPLES
	conda:
		"envs/RNASeqPipelineProject.yml"
	script:
		"scripts/make_matrix.py"

rule deseq:
	input:
		matrix="results/matrix/count_matrix.tsv",
		samples="config/samples.tsv"
	output:
		deseq="results/deseq/deseq_results.csv"
	conda:
		"envs/deseq2.yml"
	script:
		"scripts/DE.R"
		
rule volcano:
	input:
		deseq="results/deseq/deseq_results.csv"
	output:
		pdf="results/volcano/volcano.pdf",
		png="results/volcano/volcano.png",
		csv="results/volcano/deseqUpdated.csv"
	conda:
		"envs/deseq2.yml"
	script:
		"scripts/volcano.R"
		