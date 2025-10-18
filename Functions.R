needed.packages <- c("ComplexHeatmap","AnnotationDbi","circlize","corrplot","cowplot","clusterProfiler",
                     "dplyr", "DESeq2", "DOSE","data.table","EnhancedVolcano","enrichplot","forcats","ggplot2",  "gridExtra", "ggpubr", "ggnewscale","ggplotify","ggrepel", "Hmisc","janitor","lubridate","limma","purrr","PCAtools", "org.Mm.eg.db","org.Hs.eg.db","jtools",
                     "ReactomePA","RColorBrewer","stringr","scales","sva","tidyr","tibble","UpSetR","VennDiagram","gt","DEGreport","gprofiler2","readxl","readr")
for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}


expression_plot_log <- function(df,Name,binwidth){
  
  means <- aggregate(logTPM ~  Group, df, mean)
  means$logTPM <- round(means$logTPM,2)
  p <- ggplot(df, aes(x=Group, y=logTPM,fill=Group)) + theme_bw(base_size = 14) 
  p <- p + geom_dotplot(binaxis='y', stackdir='center',binwidth=binwidth) + scale_colour_identity()
  #p <- p + ggtitle(expression(bold(paste(Name))))
  p <- p + ggtitle(Name)
  p <- p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + stat_summary(fun.y=mean, geom="point",color="black", size=2)
  p <- p + ylab("Expression in log TPM") 
  p <- p + xlab("") 
  p <- p + theme(plot.title = element_text(hjust=0.5,face="bold", color="black"),axis.title.x =element_text(face="bold"),axis.title.y = element_text(size=16))
  p <- p + theme(axis.text.x= element_text(face="bold",size=12,angle=45,vjust=0.6))
  p <- p + theme(axis.text.y= element_text(size=16))
  p <- p + theme(legend.position="none")
  p <- p + geom_text(data = means, aes(label = logTPM, hjust=1.2, y = logTPM + 0.02),size=6)
  p <- p + stat_compare_means(label = "p.format",method = "t.test",comparisons = my_comparisons)
  p
}

plotPCA.df <-
  function (object,
            metadata = metadata,
            intgroup = "condition",
            ntop = 500,
            returnData = FALSE)
  {
    rv <- rowVars((object))
    select <-
      order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(object[select,]))
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    if (!all(intgroup %in% names(metadata))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(metadata[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
      metadata[[intgroup]]
    }
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        PC3 = pca$x[, 3],
        PC4 = pca$x[, 4],
        PC5 = pca$x[, 5],
        PC6 = pca$x[, 6],
        PC7 = pca$x[, 7],
        PC8 = pca$x[, 8],
        #PC9 = pca$x[, 9],
        #PC10 = pca$x[, 10] ,
        intgroup.df,
        name = rownames(metadata)
      )
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:8]
      return(d)
    }
  }

map_function.df <- function(x, inputtype, outputtype) {
  mapIds(
    org.Mm.eg.db,
    keys = row.names(x),
    column = outputtype,
    keytype = inputtype,
    multiVals = "first"
  )
}


DEG.SampleType <- function(rawdata,meta) {
  dseq_res <- data.frame()
  All_res <- data.frame()
  
  dat2 <- as.matrix(rawdata[, colnames(rawdata) %in% rownames(meta)])

  ddsHTSeq <- DESeqDataSetFromMatrix(countData = dat2,
                                     colData = meta,
                                     design = ~ SampleType)
  ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 8,]
  #ddsHTSeq$Genotype <- relevel(ddsHTSeq$Genotype, ref = "LOAD2")
  dds <- DESeq(ddsHTSeq, parallel = TRUE)
  res <- results(dds, alpha = 0.05)
  summary(res)
  res$symbol <- map_function.df(res, "ENSEMBL", "SYMBOL")
  res$EntrezGene <- map_function.df(res, "ENSEMBL", "ENTREZID")
  
  All_res <<- as.data.frame(res[, c(7:8, 1:6)])
  
}


## function to find overlap between DEG gene lists
overlap.fn <- function(gp1, DEG.genes)
{
  dat.up <- list()
  for (i in 1:length(gp1)) {
    dat.up[[i]] <-
      DEG.genes[DEG.genes$group %in% gp1[i] &
                  DEG.genes$log2FoldChange > 0,] %>% pull(symbol) %>% na.omit()
    names(dat.up)[i] <- gp1[i]
  }
  
  dat.down <- list()
  for (i in 1:length(gp1)) {
    dat.down[[i]] <-
      DEG.genes[DEG.genes$group %in% gp1[i] &
                  DEG.genes$log2FoldChange < 0,] %>% pull(symbol) %>% na.omit()
    names(dat.down)[i] <- gp1[i]
  }
  all.up <<- dat.up
  all.down <<- dat.down
}

## KEGG enrichment analysis
kegg.fn <- function(gp1, DEG.genes,n)
{
  dat.up <- list()
  for (i in 1:length(gp1)) {
    dat.up[[i]] <-
      DEG.genes[DEG.genes$group %in% gp1[i] &
                  DEG.genes$log2FoldChange > 0,] %>% pull(EntrezGene)
    names(dat.up)[i] <- gp1[i]
  }
  
  dat.down <- list()
  for (i in 1:length(gp1)) {
    dat.down[[i]] <-
      DEG.genes[DEG.genes$group %in% gp1[i] &
                  DEG.genes$log2FoldChange < 0,] %>% pull(EntrezGene)
    names(dat.down)[i] <- gp1[i]
  }
  
  ## perform enrichment analysis
  enrich_pathway.up <- compareCluster(dat.up,
                                      fun = "enrichKEGG",
                                      pvalueCutoff = 0.1,
                                      organism = "mmu")
  
  enrich_pathway.up@compareClusterResult$Description <-
    gsub(
      " - Mus musculus \\(house mouse)",
      "",
      enrich_pathway.up@compareClusterResult$Description
    )
  print(
    clusterProfiler::dotplot(
      enrich_pathway.up,
      showCategory = n,
      font.size = 12,
      label_format = 60
    ) + ggtitle("KEGG Enrichment in upregulated genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")
  )
  
  ## perform enrichment analysis
  enrich_pathway.down <- compareCluster(
    dat.down,
    fun = "enrichKEGG",
    pvalueCutoff = 0.1,
    organism = "mmu"
  )
  
  enrich_pathway.down@compareClusterResult$Description <-
    gsub(
      " - Mus musculus \\(house mouse)",
      "",
      enrich_pathway.down@compareClusterResult$Description
    )
  print(
    clusterProfiler::dotplot(
      enrich_pathway.down,
      showCategory = n,
      font.size = 12,
      label_format = 60
    ) + ggtitle("KEGG Enrichment in downregulated genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")
  )
}

# df_up <- as.data.frame(enrich_pathway.up) %>% arrange(p.adjust) %>% head(10) 
# df_up$Regulation <- "Up-regulated" 
# 
# df_down <- as.data.frame(enrich_pathway.down) %>% arrange(p.adjust) %>% head(10)
# df_down$Regulation <- "Down-regulated"
# 
# combined_df <- rbind(df_up, df_down)
# ggplot(combined_df, aes(x = -log10(p.adjust), y = Description, fill = Regulation)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "KEGG Pathway Enrichment for Up and Down Regulated Genes",
#        x = "-log10(Adjusted P-value)",
#        y = "Pathway Description") +
#   facet_grid(rows = vars(Regulation), cols=vars(Cluster),scales = "free_y", space = "free") +
#   theme_minimal()
# 
# ggplot(combined_df, aes(x = Cluster, y = Description, fill = -log10(p.adjust))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "KEGG Pathway Enrichment for Up and Down Regulated Genes",
#        x = "-log10(Adjusted P-value)",
#        y = "Pathway Description") +
#   facet_grid(~Regulation) + 
#   theme_minimal()


KEGG.GSEA <-
  function(kegg_gene_list) {
    keggpath.XX <-
      gseKEGG(
        geneList     = kegg_gene_list,
        organism     = 'mmu',
        nPerm        = 10000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "fdr",
        keyType       = "ncbi-geneid"
      )
    keggpath.XX
  }

gsea.fn <- function(dat) {
  genot3 <- unique(dat$group)
  pathways2 <- data.frame()
  for (i in 1:length(genot3))
  {
    df <- dat %>% filter(group %in% genot3[i])
    geneList <- df %>% pull(log2FoldChange)
    names(geneList) <- as.character(df %>% pull(EntrezGene))
    kegg_gene_list <- sort(geneList, decreasing = TRUE)
    kk2 <- KEGG.GSEA(kegg_gene_list)
    kk3 <- cbind(as.data.frame(kk2), Genotype = genot3[i])
    pathways2 <- rbind(pathways2, kk3)
  }
  
  GSEA.PS <- pathways2
  GSEA.PS$Description <-
    gsub(" - Mus musculus \\(house mouse)", "", GSEA.PS$Description)
  GSEA.PS$Genotype <- factor(GSEA.PS$Genotype, levels = genot3)
  return(GSEA.PS)
}


magora_corrplot_colored <- function(data, ran)
{
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$module, y = .data$model_sex),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$module,
        y = .data$model_sex,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      ggplot2::aes(
        x = .data$module,
        y = .data$model_sex,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Model(Sex)| Age/Control") +
    ggplot2::facet_grid(
      #rows = dplyr::vars(.data$group),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = 'y'
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 14, colour = c("black")),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 14),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 14
      ),
      axis.text.y = ggplot2::element_text(angle = 0, size = 14),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -54,
        hjust = 0.03,
        size = 12,
        face = "bold"
      ),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) + theme(legend.title = element_text(size=13),legend.text = element_text(size=12))
  
  fills <- c("darkorange3","chartreuse3","deepskyblue2","turquoise","deeppink2")
  g <- ggplot_gtable(ggplot_build(p2))
  stripr <- which(grepl('strip-t', g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  ggdraw(g)
}

magora_corrplot <- function(data, ran)
{
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$module, y = .data$model_sex),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$module,
        y = .data$model_sex,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      ggplot2::aes(
        x = .data$module,
        y = .data$model_sex,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Model(Sex)| Age/Control") +
    ggplot2::facet_grid(
      rows = dplyr::vars(.data$group),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = 'y'
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 14, colour = c("black")),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 14),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 14
      ),
      axis.text.y = ggplot2::element_text(angle = 0, size = 14),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -54,
        hjust = 0.03,
        size = 12,
        face = "bold"
      ),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

corr_function_lm <- function(data,human_data)
{
  ns_vs_ampad_fc <- data  %>%
    inner_join(human_data, by = c("Gene")) %>%
    group_by(module, Variant) %>%
    nest(data = c(Gene, value, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["value"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  
  # Process data for plotting ----
  # Flag for significant results, add cluster information to modules
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < 0.05, age_group = "All Months") %>%
    left_join(module_clusters, by = "module") %>%
    dplyr::select(
      cluster,
      cluster_label,
      module,
      Variant,
      age_group,
      correlation = estimate,
      p_value,
      significant
    )
  
  # Create a version of the data for plotting - clean up naming, order factors, etc
  nanostring_for_plot.all <- nanostring %>%
    arrange(cluster) %>%
    mutate(
      Variant = factor(Variant, levels = ordered.variant),
      Variant = fct_rev(Variant),
      module = factor(module, levels = mod),
    )
}

variant_corrplot <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$module, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::facet_grid(
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 9),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(face = "italic", size = 12),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

variant_corrplot_test <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$module, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(
      rows = dplyr::vars(.data$Background),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 11),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 12),
      strip.background.y = ggplot2::element_rect(fill = "grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -56,
        hjust = 0.01,
        size = 11,
        face = "bold"
      ),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

variant_corrplot_colored <- function(data, ran) {
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$module, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$module,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(
      rows = dplyr::vars(.data$Background),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 11),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 12),
      strip.background.y = ggplot2::element_rect(fill = "grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -56,
        hjust = 0.01,
        size = 11,
        face = "bold"
      ),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) + theme(legend.title = element_text(size=13),legend.text = element_text(size=12))
  
  fills <- c("darkorange3","chartreuse3","deepskyblue2","turquoise","deeppink2")
  g <- ggplot_gtable(ggplot_build(p2))
  stripr <- which(grepl('strip-t', g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  ggdraw(g)
  
}


subtypecorr_function_lm <- function(data)
{
  ns_vs_ampad_fc <- data  %>%
    inner_join(ampad_subtype_fc, by = c("Gene")) %>%
    group_by(subtype, Variant) %>%
    nest(data = c(Gene, value, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["value"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  
  # Process data for plotting ----
  # Flag for significant results, add cluster information to modules
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < 0.05, age_group = "All months") %>%
    left_join(module_cohort, by = "subtype") %>%
    dplyr::select(
      cluster,
      cluster_label,
      subtype,
      Variant,
      age_group,
      correlation = estimate,
      p_value,
      significant
    )
  
  # Create a version of the data for plotting - clean up naming, order factors, etc
  nanostring_for_plot.all <- nanostring %>%
    arrange(cluster) %>%
    mutate(
      Variant = factor(Variant, levels = ordered.variant),
      Variant = fct_rev(Variant),
      subtype = factor(subtype, levels = subtypes),
    )
}

Neffsubtypecorr_function_lm <- function(data)
{
  ns_vs_ampad_fc <- data  %>%
    inner_join(neff_subtype_fc, by = c("Gene")) %>%
    group_by(subtype, Variant) %>%
    nest(data = c(Gene, value, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["value"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  
  # Process data for plotting ----
  # Flag for significant results, add cluster information to modules
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < 0.05, age_group = "All months") %>%
    #left_join(module_cohort, by = "subtype") %>%
    dplyr::select(subtype,
                  Variant,
                  age_group,
                  correlation = estimate,
                  p_value,
                  significant) %>% mutate(cluster_label = "MSBB_PHG_Neff")
  # Create a version of the data for plotting - clean up naming, order factors, etc
  nanostring_for_plot.all <- nanostring %>%
    #arrange(cluster) %>%
    mutate(
      Variant = factor(Variant, levels = ordered.variant),
      Variant = fct_rev(Variant)
      #subtype = factor(subtype,levels=subtypes),
    )
}

subtype_variant_corrplot <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$subtype, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::facet_grid(
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 9),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(face = "italic", size = 12),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

subtype_variant_corrplot_test <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$subtype, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(
      rows = dplyr::vars(.data$age_group),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 12),
      strip.background.y = ggplot2::element_rect(fill = "grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(face = "italic", size = 12),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -56,
        hjust = 0.04,
        size = 11,
        face = "bold"
      ),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

subtype_variant_corrplot_test2 <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$subtype, y = .data$Variant),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      aes(
        x = .data$subtype,
        y = .data$Variant,
        colour = .data$correlation
      ),
      color = "black",
      shape = 0,
      size = 9
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran),
      breaks = c(-ran, 0, ran),
      low = "#85070C",
      high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(
      rows = dplyr::vars(.data$age_group),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 12),
      strip.background.y = ggplot2::element_rect(fill = "grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = 12
      ),
      axis.text.y = ggplot2::element_text(face = "italic", size = 12),
      plot.title = ggplot2::element_text(
        angle = 0,
        vjust = -40,
        hjust = 0.04,
        size = 11,
        face = "bold"
      ),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

corr_function_lm_biodomain <- function(data,ref,sig,age)
{
  ns_vs_ampad_fc <- data  %>% 
    inner_join(ref, by = c("Gene")) %>%
    group_by(Term,Variant) %>%
    nest(data = c(Gene, value, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["value"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  
  # Process data for plotting ----
  # Flag for significant results, add cluster information to modules
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < sig,age_group=age) %>%
    dplyr::select(biodom=Term, Variant,age_group, correlation = estimate, p_value, significant)
  
  # Create a version of the data for plotting - clean up naming, order factors, etc
  nanostring_for_plot.all <- nanostring %>%
    arrange(biodom) %>%
    mutate(
      Variant =factor(Variant,levels=ordered.variant),
      Variant =fct_rev(Variant),
      module = factor(biodom,levels=domain),
    ) 
}

variant_corrplot_biodom <- function(data,ran) {
  
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = data, ggplot2::aes(x = .data$biodom, y = .data$Variant), colour = "black", fill = "white") +
    ggplot2::geom_point(data = dplyr::filter(data), ggplot2::aes(x = .data$biodom, y = .data$Variant, colour = .data$correlation, size = abs(.data$correlation))) +
    ggplot2::geom_point(data = dplyr::filter(data, .data$significant),aes(x=.data$biodom,y=.data$Variant, colour = .data$correlation),color="black",shape=0,size=9) +
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) + 
    ggplot2::scale_color_gradient2(limits = c(-ran, ran), breaks = c(-ran, 0, ran), low = "#85070C", high = "#164B6E", name = "Correlation", guide = ggplot2::guide_colorbar(ticks = FALSE)) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(rows = dplyr::vars(.data$group),cols = dplyr::vars(.data$brain),scales = "free", space = "free",switch="y") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(angle = 0,size = 14),
      strip.text.y.left = ggplot2::element_text(angle = 0,size = 14),
      strip.background.y = ggplot2::element_rect(fill="grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0,size=14),
      axis.text.y = ggplot2::element_text(size=14),
      plot.title = ggplot2::element_text(angle = 0, vjust = -56, hjust = 0.01,size=11,face="bold"),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

human_biodomains <- tibble(
  Biodomain = c("Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair", 
                "Endolysosome", "Epigenetic", "Immune Response", "Lipid Metabolism", 
                "Metal Binding and Homeostasis", "Mitochondrial Metabolism", 
                "Myelination", "Oxidative Stress", "Proteostasis", "RNA Spliceosome", 
                "Structural Stabilization", "Synapse", "Tau Homeostasis", "Vasculature"),
  Biodomain_label = c("Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair", 
                      "Endolysosome", "Epigenetic", "Immune\n Response", "Lipid Metabolism", 
                      "Metal Binding\n and Homeostasis", "Mitochondrial\n Metabolism", 
                      "Myelination", "Oxidative Stress", "Proteostasis", "RNA\n Spliceosome", 
                      "Structural Stabilization", "Synapse", "Tau\n Homeostasis", "Vasculature")
)

corr_function_lm_subbiodomain <- function(data,ref,sig,age)
{
  ns_vs_ampad_fc <- data  %>% 
    inner_join(ref, by = c("Gene")) %>%
    group_by(Term,Variant) %>%
    nest(data = c(Gene, value, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["value"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  
  # Process data for plotting ----
  # Flag for significant results, add cluster information to modules
  nanostring <- ns_vs_ampad_fc %>%
    left_join(human_biodomains, by = "Biodomain") %>%
    mutate(significant = p_value < 0.05,age_group=age) %>%
    dplyr::select(module=Biodomain,biodom=Term, Variant,age_group, correlation = estimate, p_value, significant,Biodomain_label)
  
  nanostring_for_plot.all <- nanostring %>%
    arrange(module) %>%
    mutate(
      Variant =factor(Variant,levels=ordered.variant),
      Variant =fct_rev(Variant),
      module = factor(module,levels=domain)
    )
}

variant_corrplot_subbiodom <- function(data,ran) {
  
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = data, ggplot2::aes(x = .data$biodom, y = .data$Variant), colour = "black", fill = "white") +
    ggplot2::geom_point(data = dplyr::filter(data), ggplot2::aes(x = .data$biodom, y = .data$Variant, colour = .data$correlation, size = abs(.data$correlation))) +
    ggplot2::geom_point(data = dplyr::filter(data, .data$significant),aes(x=.data$biodom,y=.data$Variant, colour = .data$correlation),color="black",shape=0,size=9) +
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) + 
    ggplot2::scale_color_gradient2(limits = c(-ran, ran), breaks = c(-ran, 0, ran), low = "#85070C", high = "#164B6E", name = "Correlation", guide = ggplot2::guide_colorbar(ticks = FALSE)) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(rows = dplyr::vars(.data$group),cols = dplyr::vars(.data$Biodomain_label), scales="free", space = "free",switch="y") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(angle = 0,size = 14),
      strip.text.y.left = ggplot2::element_text(angle = 0,size = 14),
      strip.background.y = ggplot2::element_rect(fill="grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0,size=14),
      axis.text.y = ggplot2::element_text(size=12),
      plot.title = ggplot2::element_text(angle = 0, vjust = -56, hjust = 0.01,size=11,face="bold"),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

variant_corrplot_subbiodom.Potrait <- function(data,ran) {
  
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = data, ggplot2::aes(y = .data$biodom, x = .data$Variant), colour = "black", fill = "white") +
    ggplot2::geom_point(data = dplyr::filter(data), ggplot2::aes(y = .data$biodom, x = .data$Variant, colour = .data$correlation, size = abs(.data$correlation))) +
    ggplot2::geom_point(data = dplyr::filter(data, .data$significant),aes(y=.data$biodom,x=.data$Variant, colour = .data$correlation),color="black",shape=0,size=9) +
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::scale_size(guide = "none", limits = c(0, ran)) + 
    ggplot2::scale_color_gradient2(limits = c(-ran, ran), breaks = c(-ran, 0, ran), low = "#85070C", high = "#164B6E", name = "Correlation", guide = ggplot2::guide_colorbar(ticks = FALSE)) +
    ggplot2::labs(x = NULL, y = NULL) +
    #ggplot2::ggtitle("Perturbation| Control") +
    ggplot2::facet_grid(cols = dplyr::vars(.data$group),rows = dplyr::vars(.data$Biodomain_label), scales="free", space = "free",switch="y") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(angle = 0,size = 14),
      strip.text.y.left = ggplot2::element_text(angle = 0,size = 14),
      strip.background.y = ggplot2::element_rect(fill="grey95"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0,size=14),
      axis.text.y = ggplot2::element_text(size=12),
      plot.title = ggplot2::element_text(angle = 0, vjust = -56, hjust = 0.01,size=11,face="bold"),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}




enrichORA_BD <- function(dat,set1){
  
  for(i in 1:length(set1)){
    print(set1[i])
    gene.list.up <- dat %>%
      filter(group == set1[i],
             log2FoldChange > 0) %>%
      arrange(desc(log2FoldChange)) %>%
      filter(!duplicated(symbol),!is.na(symbol)) %>%
      pull(symbol) %>%
      unique()
    
    gene.list.dn <- dat %>%
      filter(group == set1[i],
             log2FoldChange < 0) %>%
      arrange(desc(log2FoldChange)) %>%
      filter(!duplicated(symbol),!is.na(symbol)) %>%
      pull(symbol) %>%
      unique()
    
    enr.up <- enrichGO(gene.list.up, 
                       ont = 'all', 
                       OrgDb = org.Mm.eg.db, 
                       keyType = 'SYMBOL',
                       universe = univ
    )
    enr.dn <- enrichGO(gene.list.dn, 
                       ont = 'all', 
                       OrgDb = org.Mm.eg.db, 
                       keyType = 'SYMBOL',
                       universe = univ
    )
    
    if(nrow(enr.dn@result) == 0 & nrow(enr.up@result) == 0){
      message("No enrichment found.")
    }
    else{
      enr.ora = bind_rows(enr.up@result %>% mutate(dir = 'up'),
                          enr.dn@result %>% mutate(dir = 'dn')) %>%
        left_join(., biodom_mouse %>% dplyr::select(Biodomain, Subdomain,ID = GO_ID), by = 'ID') %>%
        relocate(Biodomain)
      
      enr.ora$Biodomain[is.na(enr.ora$Biodomain)] <- 'none'
      enr.ora <- enr.ora %>% mutate(Subdomain = ifelse(Subdomain %in% NA,.$Biodomain,.$Subdomain))
      
      bd.tally = tibble(domain = c(unique(biodom_mouse$Biodomain), 'none')) %>%
        rowwise() %>%
        mutate(
          n_term = biodom_mouse$GO_ID[biodom_mouse$Biodomain == domain] %>% unique() %>% length(),
          n_sig_term = enr.ora$ID[enr.ora$Biodomain == domain] %>% unique() %>% length()
        )
      
      enr.ora <-
        left_join(enr.ora, dom.lab, by = c('Biodomain' = 'domain')) %>%
        mutate(Biodomain = factor(Biodomain, levels = arrange(bd.tally, n_sig_term) %>% pull(domain))) %>% mutate(model=set1[i]) %>%
        arrange(Biodomain, p.adjust) %>%
        mutate(
          signed_logP = -log10(p.adjust),
          signed_logP = if_else(dir == 'dn',-1 * signed_logP, signed_logP)
        )
      
      enr.ora.all <- bind_rows(enr.ora.all,enr.ora)
    }
    
  }
  return(enr.ora.all)
}


enrichGSEA_BD <- function(dat,set1)
{
  for(i in 1:length(set1))
  {
    gene.list <- dat %>% 
      filter(group == set1[i]) %>%
      arrange(desc(log2FoldChange)) %>% 
      filter(!duplicated(symbol), !is.na(symbol)) %>% 
      pull(log2FoldChange, name = symbol)
    
    enr <- gseGO(gene.list, ont = 'all', OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',pvalueCutoff = 1,seed=TRUE)
    
    enr.bd = enr@result %>% 
      left_join(., biodom_mouse %>% dplyr::select(Biodomain, Subdomain,ID = GO_ID), by = 'ID') %>% 
      relocate(Biodomain)
    
    enr.bd$Biodomain[is.na(enr.bd$Biodomain)] <- 'none'
    enr.bd <- enr.bd %>% mutate(Subdomain = ifelse(Subdomain %in% NA,.$Biodomain,.$Subdomain))
    
    bd.tally = tibble(domain = c(unique(biodom_mouse$Biodomain), 'none')) %>% 
      rowwise() %>% 
      mutate(
        n_term = biodom_mouse$GO_ID[ biodom_mouse$Biodomain == domain ] %>% unique() %>% length(),
        n_sig_term = enr.bd$ID[ enr.bd$Biodomain == domain ] %>% unique() %>% length()
      )
    #arrange(bd.tally, desc(n_sig_term))
    
    enr.bd <- full_join(enr.bd, dom.lab, by = c('Biodomain' = 'domain')) %>% 
      mutate(Biodomain = factor(Biodomain, levels = arrange(bd.tally, n_sig_term) %>% pull(domain))) %>% mutate(model=set1[i]) %>%
      arrange(Biodomain, p.adjust) 
    
    enr.GSEA.bd.all <- bind_rows(enr.GSEA.bd.all,enr.bd)
  }
  return(enr.GSEA.bd.all)
}

selected.path2 <- c("Osteoclast differentiation","Lysosome","Phagosome","ECM-receptor interaction","Neuroactive ligand-receptor interaction",
                    "Fc gamma R-mediated phagocytosis","Axon guidance","Focal adhesion","Cholinergic synapse","GABAergic synapse",
                    "Oxytocin signaling pathway","Wnt signaling pathway","Endocytosis","Alzheimer disease","Oxidative phosphorylation")


selected.path <- c("Alzheimer disease", "Phagosome", "Wnt signaling pathway", "Oxytocin signaling pathway", "Chagas disease", 
                   "C-type lectin receptor signaling pathway", "Dopaminergic synapse", "Proteasome",
                   "Neurotrophin signaling pathway", "cAMP signaling pathway", "Lysosome", "Axon guidance", "Ras signaling pathway", "Synaptic vesicle cycle", 
                   "Estrogen signaling pathway", "Hippo signaling pathway", "Cellular senescence", "Regulation of actin cytoskeleton", 
                   "Rap1 signaling pathway", "HIF-1 signaling pathway", "Calcium signaling pathway", "Non-alcoholic fatty liver disease", 
                   "Cholinergic synapse", "Osteoclast differentiation", "Fc gamma R-mediated phagocytosis", "Adherens junction",
                   "MAPK signaling pathway", "Natural killer cell mediated cytotoxicity", 
                   "Chemokine signaling pathway", "Glutamatergic synapse", "Thyroid hormone signaling pathway", 
                   "Insulin resistance", "Kaposi sarcoma-associated herpesvirus infection", 
                   "ECM-receptor interaction", "Huntington disease", "Apelin signaling pathway", 
                   "Proteoglycans in cancer", "Insulin signaling pathway", "Platelet activation", "cGMP-PKG signaling pathway", 
                   "Thermogenesis", "Toxoplasmosis", "Tight junction", "Oxidative phosphorylation", 
                   "AGE-RAGE signaling pathway in diabetic complications", "Epstein-Barr virus infection", 
                   "Parkinson disease", "Human cytomegalovirus infection", "Th17 cell differentiation", 
                   "Cytokine-cytokine receptor interaction", "Cell adhesion molecules", "Protein processing in endoplasmic reticulum", 
                   "Vascular smooth muscle contraction", "Focal adhesion", "Autophagy - animal", "EGFR tyrosine kinase inhibitor resistance", 
                   "Transcriptional misregulation in cancer", "Leukocyte transendothelial migration", "Human T-cell leukemia virus 1 infection", 
                   "Retrograde endocannabinoid signaling", "Phosphatidylinositol signaling system", 
                   "GABAergic synapse", "Pathways in cancer", "Endocytosis", "Adrenergic signaling in cardiomyocytes", 
                   "FoxO signaling pathway", "Morphine addiction", "Sphingolipid signaling pathway", 
                   "Neuroactive ligand-receptor interaction", "Human papillomavirus infection", 
                   "Fluid shear stress and atherosclerosis", "Metabolic pathways", "PI3K-Akt signaling pathway")
