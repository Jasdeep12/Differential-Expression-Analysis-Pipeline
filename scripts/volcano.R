suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ggplot2))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))


de_df <- read.csv(snakemake@input[["deseq"]])

de_df <- de_df[!is.na(de_df$padj), ]

de_df$symbol <- mapIds(
	org.Hs.eg.db,
	keys = de_df$X,
	keytype = "ENSEMBL",
	column = "SYMBOL",
	multiVals = "first"
	)
write.csv(de_df, snakemake@output[["csv"]])

sig <- de_df[!is.na(de_df$symbol) &
             de_df$padj < 0.05 &
             abs(de_df$log2FoldChange) > 1.0, ]
top_genes <- sig$symbol[order(sig$padj)][1:15]	
	

p <- EnhancedVolcano(de_df,
	lab = de_df$symbol,
	selectLab = top_genes,
	x = 'log2FoldChange',
	y = 'padj',
	pCutoff = 0.05,
	FCcutoff = 1.0,
	pointSize = 2.0,
	labSize = 4.0,
	drawConnectors = TRUE,
	title = 'Human bronchial epithielial cells: Normal VS Cystic Fibrosis',
	subtitle = 'Differential expression analysis',
	legendPosition = 'right')

ggsave(snakemake@output[["pdf"]], plot = p, width = 14, height = 8, units = "in", dpi = 300)

ggsave(snakemake@output[["png"]], plot = p, width = 14, height = 8, units = "in", dpi = 300)






