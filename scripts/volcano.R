suppressMessages(library(EnhancedVolcano))
suppressMessagesq(library(ggplot2))


de <- read.delim(snakemake@input[["deseq"]])
de_df <- as.data.frame(de)

de_df <- de_df[!is.na(de_df$padj), ]
