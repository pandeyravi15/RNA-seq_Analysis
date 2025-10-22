### KEGG Pathways Enrichment

#Over-representation (or enrichment) analysis is a statistical method that determines whether genes from pre-defined sets (ex: those beloging to a specific GO term or KEGG pathway) are present more than would be expected (over-represented) in a subset of your data. In this case, the subset is your set of significantly under or over expressed genes identified in DESeq2 analysis. 
#We look for enrichment of biological pathways in a list of differentially expressed genes. Here we test for enrichment of KEGG pathways using using enrichKEGG function in [clusterProfiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html) package and plotted enriched terms at `p.adjust < 0.1`. 

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

dat <- DE_Genotype.df 
DEG.genes <- subset(dat[order(dat$padj),],padj<0.05)
gp1 <- unique((DEG.genes$group))
kegg.fn(gp1,DEG.genes,n=15)

## GSEA Analysis
dat <- DE_Genotype.df %>% select(EntrezGene, log2FoldChange, group) %>% na.omit()

genot3 <- unique(dat$group)
pathways2 <- data.frame()
for (i in 1:length(genot3))
{
  df <- dat %>% filter(group %in% genot3[i]) %>% distinct(.,EntrezGene,.keep_all = TRUE)
  geneList <- df %>% pull(log2FoldChange)
  names(geneList) <- df %>% pull(EntrezGene)
  kegg_gene_list <- sort(geneList, decreasing = TRUE)
  kk2 <- KEGG.GSEA(kegg_gene_list)
  kk3 <- cbind(as.data.frame(kk2), Genotype = genot3[i])
  pathways2 <- rbind(pathways2, kk3)
}

GSEA.PS <- pathways2
GSEA.PS$Description <-
  gsub(" - Mus musculus \\(house mouse)", "", GSEA.PS$Description)
GSEA.PS$Genotype <- factor(GSEA.PS$Genotype, levels = genot3)

# save results
#save(GSEA.PS,file="../results/GSEA_LOAD3_Brain_Transcriptomics.RData")

# load("../results/GSEA_LOAD3_Brain_Transcriptomics.RData")
mat1 <-
  GSEA.PS[(GSEA.PS$Description %in% selected.path), ] %>% dplyr::select(Genotype, Description, NES) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Genotype",
              values_from = "NES") %>% column_to_rownames(var = "Description")  %>% as.data.frame(.) %>% na.omit()

mat2 <-
  GSEA.PS[(GSEA.PS$Description %in% selected.path), ] %>% dplyr::select(Genotype, Description, p.adjust) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Genotype",
              values_from = "p.adjust") %>% column_to_rownames(var = "Description")  %>% as.data.frame(.) %>% na.omit()

mat3 <-
  GSEA.PS[(GSEA.PS$Description %in% selected.path), ] %>% dplyr::select(Genotype, Description, pvalue) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Genotype",
              values_from = "pvalue") %>% column_to_rownames(var = "Description")  %>% as.data.frame(.) %>% na.omit()

## check term names ar ein order in both matrix or not
#match(rownames(mat1),rownames(mat2))
#all(rownames(mat1) == rownames(mat2))
#all(rownames(mat1) == rownames(mat3))

col_fun = colorRamp2(c(-2, 0, 2), c("red", "white", "blue"))
Heatmap(
  as.matrix((mat1)),
  name = "NES",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_row_names = TRUE,
  row_names_max_width = unit(11, "cm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 18, fontface = "italic"),
  row_names_gp = gpar(fontsize = 14),
  column_names_rot = 60,
  column_names_max_height = unit(12, "cm"),
  col = col_fun,
  heatmap_legend_param = list(
    legend_height = unit(6, "cm"),
    grid_width = unit(1, "cm"),
    title_gp = gpar(fontsize = 16),
    labels_gp = gpar(fontsize = 16)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(mat2[i, j] < 0.1 & mat2[i, j] > 0.05) {
      grid.text("*", x, y)
    } else if(mat2[i, j] < 0.05) {
      grid.text("**", x, y)
    }
    
  }
)

Heatmap(
  as.matrix((mat1)),
  name = "NES",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  show_row_names = TRUE,
  row_names_max_width = unit(11, "cm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 18, fontface = "italic"),
  row_names_gp = gpar(fontsize = 14),
  column_names_rot = 60,
  column_names_max_height = unit(12, "cm"),
  col = col_fun,
  heatmap_legend_param = list(
    legend_height = unit(6, "cm"),
    grid_width = unit(1, "cm"),
    title_gp = gpar(fontsize = 16),
    labels_gp = gpar(fontsize = 16)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(mat3[i, j] < 0.1 & mat2[i, j] > 0.05) {
      grid.text("*", x, y)
    } else if(mat3[i, j] < 0.05) {
      grid.text("**", x, y)
    }
    
  }
)