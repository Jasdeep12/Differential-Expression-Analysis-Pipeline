suppressMessages({Library("clusterProfiler"),
Library("org.Hs.eg.db"),
Library("enrichplot"),
Library("ggplot2"),
Library("msigdbr"),
})

de <- read.csv(snakemake@param[["de"]])

set.seed(1234)

ranked <- de %>%
    filter(!is.na(padj)) %>%
    arrange(desc(stat)) %>%
    pull(stat, name = G)

mark <- msigdbr(collection = "H") %>%
    dplyr::select(gs_name, ensembl_gene)

mark2 <- msigdbr(collection = "H") %>%
    dplyr::select(gs_name, gs_description)

eh <- GSEA(ranked, TERM2GENE = mark, TERM2NAME = mark2, seed = TRUE)

