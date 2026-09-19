#scripts/DE.R

suppressMessages(library(DESeq2))


# Loading Inputs
cts <- read.delim(snakemake@input[["matrix"]], row.names = "Geneid")
sinfo <- read.delim(snakemake@input[["samples"]])


# Align Sample Order
# Sample order comes pre aligned if system works properly
stopifnot(all(colnames(cts) == sinfo$sample))

sinfo$condition <- relevel(factor(sinfo$condition), ref = "control")


# Building DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
	countData = cts,
	colData = sinfo,
	design = ~ condition
	)

# Running DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Write Results
write.csv(as.data.frame(res),snakemake@output[[1]])





