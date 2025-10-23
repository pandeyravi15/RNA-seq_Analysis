### KEGG Pathways Enrichment
#Over-representation (or enrichment) analysis is a statistical method that determines whether genes from pre-defined sets (ex: those beloging to a specific GO term or KEGG pathway) are present more than would be expected (over-represented) in a subset of your data. In this case, the subset is your set of significantly under or over expressed genes identified in DESeq2 analysis. 
#We look for enrichment of biological pathways in a list of differentially expressed genes. Here we test for enrichment of KEGG pathways using using enrichKEGG function in [clusterProfiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html) package and plotted enriched terms. 

# LOAD the results from DESeq analysis
lnames  = load("data/DESeq_Results_Transcripotmics.RData")
DEG.genes <- subset(DE_Genotype.df[order(DE_Genotype.df$padj),],padj<0.05)
gp1 <- unique((DEG.genes$group))

# storing all up-regulated genes in a list for each comparison
dat.up <- list()
for (i in 1:length(gp1)) {
  dat.up[[i]] <-
    DEG.genes[DEG.genes$group %in% gp1[i] &
                DEG.genes$log2FoldChange > 0,] %>% pull(EntrezGene)
  names(dat.up)[i] <- gp1[i]
}
 
# storing all down-regulated genes in a list for each comparison
dat.down <- list()
for (i in 1:length(gp1)) {
  dat.down[[i]] <-
    DEG.genes[DEG.genes$group %in% gp1[i] &
                DEG.genes$log2FoldChange < 0,] %>% pull(EntrezGene)
  names(dat.down)[i] <- gp1[i]
}
  
## perform enrichment analysis for KEGG pathways
enrich_pathway.up <- compareCluster(
  dat.up,
  organism = "mmu",
  fun = "enrichKEGG",
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH"
)


enrich_pathway.down <- compareCluster(
  dat.down,
  organism = "mmu",
  fun = "enrichKEGG",
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH"
)


## plot top enriched KEGG pathways
print(
  clusterProfiler::dotplot(
    enrich_pathway.up,
    showCategory = 10,
    font.size = 12,
    label_format = 60
  ) + 
    ggtitle("KEGG Enrichment in upregulated genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")
)

print(
  clusterProfiler::dotplot(
    enrich_pathway.down,
    showCategory = 10,
    font.size = 12,
    label_format = 60
  ) + 
    ggtitle("KEGG Enrichment in downregulated genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")
)

## Gene Set Enrichment Analysis
dat <- DE_Genotype.df %>% select(EntrezGene, log2FoldChange, group) %>% na.omit()
groups <- unique(dat$group)
gsea_pathways <- data.frame()

for (i in 1:length(groups))
{
  df <- dat %>% filter(group %in% groups[i]) %>% distinct(.,EntrezGene,.keep_all = TRUE)
  geneList <- df %>% pull(log2FoldChange)
  names(geneList) <- df %>% pull(EntrezGene)
  kegg_gene_list <- sort(geneList, decreasing = TRUE)
  gsea_result <- gseKEGG(
    geneList     = kegg_gene_list,
    organism     = 'mmu',
    nPerm        = 10000,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "fdr",
    seed = TRUE,
    keyType       = "ncbi-geneid"
  )
  df <- cbind(as.data.frame(gsea_result), Genotype = groups[i])
  gsea_pathways <- rbind(gsea_pathways, df)
}

gsea_pathways <- gsea_pathways %>% mutate(Genotype = factor(Genotype, levels = groups))

# save results
#save(gsea_pathways,file="../results/GSEA_Transcriptomics.RData")

# Plot the results for selected pathways
selected.path <- c("Osteoclast differentiation","Lysosome","Phagosome","ECM-receptor interaction","Neuroactive ligand-receptor interaction",
                   "Fc gamma R-mediated phagocytosis","Axon guidance","Focal adhesion","Cholinergic synapse","GABAergic synapse",
                   "Oxytocin signaling pathway","Wnt signaling pathway","Endocytosis","Alzheimer disease","Oxidative phosphorylation")

mat1 <-
  gsea_pathways[(gsea_pathways$Description %in% selected.path), ] %>% dplyr::select(Genotype, Description, NES) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Genotype",
              values_from = "NES") %>% column_to_rownames(var = "Description")  %>% as.data.frame(.) %>% na.omit()

mat2 <-
  gsea_pathways[(gsea_pathways$Description %in% selected.path), ] %>% dplyr::select(Genotype, Description, p.adjust) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Genotype",
              values_from = "p.adjust") %>% column_to_rownames(var = "Description")  %>% as.data.frame(.) %>% na.omit()


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

## GO-terms enrichment
univ <- as.data.frame(org.Mm.egGO) %>% 
  pull(gene_id) %>% 
  unique() %>% 
  bitr(., fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = org.Mm.eg.db, drop = T) %>% 
  pull('SYMBOL') %>% 
  intersect(., DE_Genotype.df$symbol)

#Now let's test for enriched GO terms (this can take 3-4 minutes)
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