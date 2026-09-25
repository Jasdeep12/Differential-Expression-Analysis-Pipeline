suppressMessages({
library("clusterProfiler")
library("org.Hs.eg.db")
library("tidyverse")
library("ggplot2")
library("msigdbr")
})

de <- read.csv(snakemake@input[["deseq"]])
species <- snakemake@config[["species"]][[1]]
metric <- snakemake@config[["gseaStat"]][[1]]
geneSet <- snakemake@config[["geneSet"]][[1]]

set.seed(1234)

de$X <- sub("\\.\\d+$", "", de$X)

ranked <- de %>%
    filter(!is.na(.data[[metric]])) %>%
    arrange(desc(.data[[metric]])) %>%
    pull(.data[[metric]], name = X)

mark <- msigdbr(species = species, collection = geneSet) %>%
    dplyr::select(gs_name, ensembl_gene)

gsea_result <- GSEA(ranked, TERM2GENE = mark,  seed = FALSE)

write.csv(as.data.frame(gsea_result), snakemake@output[["gsea"]], row.names = FALSE)


gsea_df <- as.data.frame(gsea_result) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    mutate(
        Description = factor(Description, levels = rev(Description)),
        direction = ifelse(NES > 0, "Enriched in CF", "Enriched in control")
    )

plot_title <- snakemake@config[["plot"]][["barplot"]]
p <- ggplot(gsea_df, aes(x = NES, y = Description, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("Enriched in CF" = "lightcoral", "Enriched in control" = "lightblue")) +
    labs(
        title = plot_title,
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        fill = NULL
    ) +
    theme_minimal(base_size = 12)

ggsave(snakemake@output[["barplot"]], plot = p, width = 10, height = 8, units = "in", dpi = 300)



