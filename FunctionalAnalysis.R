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
