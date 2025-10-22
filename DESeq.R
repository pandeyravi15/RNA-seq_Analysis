## function for DESeq analysis
DEG <- function(countdata,
                meta,
                include.batch = FALSE,
                ref = ref) {
  dseq_res <- data.frame()
  All_res <- data.frame()
  
  if (include.batch) {
    cat("Including batch as covariate\n")
    design_formula <- ~ Batch + SampleType
  }
  else{
    design_formula <- ~ SampleType
  }
  
  dat2 <- as.matrix(rawdata[, colnames(rawdata) %in% rownames(meta)])
  ddsHTSeq <- DESeqDataSetFromMatrix(countData = dat2,
                                     colData = meta,
                                     design = ~ design_formula)
  ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 6,]
  ddsHTSeq$SampleType <- relevel(ddsHTSeq$SampleType, ref = ref)
  dds <- DESeq(ddsHTSeq, parallel = TRUE)
  res <- results(dds, alpha = 0.05)
  summary(res)
  res$symbol <- map_function.df(res, "ENSEMBL", "SYMBOL")
  res$EntrezGene <- map_function.df(res, "ENSEMBL", "ENTREZID")
  
  All_res <<- as.data.frame(res[, c(7:8, 1:6)])
  
}


comparisons <-  data.frame(control=c("C57BL/6J-Female-4M" ,  "C57BL/6J-Male-4M" , "C57BL/6J-Female-12M"  ,  "C57BL/6J-Male-12M"),
                           case=c("B6J.LOAD3-Female-4M", "B6J.LOAD3-Male-4M",   "B6J.LOAD3-Female-12M", "B6J.LOAD3-Male-12M" ))
# # genotype + diet
DE_Genotype.list <- list()
DE_Genotype.df <- data.frame()
#
for (i in 1:nrow(comparisons)){
  meta <- metadata[metadata$Group %in% comparisons[i,] ,] %>% rename("SampleType" = "Genotype")
  DEG.SampleType(rawcountdata,meta)
  
  #append results in data frame
  DE_Genotype.df <- rbind(DE_Genotype.df,All_res %>% mutate(model=gsub("-.*$","",comparisons[i,2])[1],sex=unique(meta$Sex),Age=sapply(strsplit(comparisons[i,2], "-"), "[", 3), group = paste0( model,"-B6","(",sex,"-",Age,")"))
  )
  
  #append results in list
  DE_Genotype.list[[i]] <- All_res
  names(DE_Genotype.list)[i] <- paste0(model=gsub("-.*$","",comparisons[i,2])[1],"-", "B6","(",sex=unique(meta$Sex),"-",Age=sapply(strsplit(comparisons[i,2], "-"), "[", 3),")")
}

#save(DE_Genotype.list,DE_Genotype.df,file="~/results/DESeq_Results_Transcripotmics.RData")


## Summary table of differential analysis results
degs.up <- map(DE_Genotype.list, ~length(which(.x$padj<0.05 & .x$log2FoldChange>0)))
degs.down <- map(DE_Genotype.list, ~length(which(.x$padj<0.05 & .x$log2FoldChange<0)))
deg1 <- data.frame(comparison=names(degs.up), Up_DEGs.pval.05=as.vector(unlist(degs.up)),Down_DEGs.pval.05=as.vector(unlist(degs.down)))

deg1 %>% gt() %>%
  tab_style(
    style = list(cell_fill(color = "lightblue"),
                 cell_text(style = "italic")),
    locations = cells_body(rows = comparison %like% "Female")
  ) %>% tab_options(table.width = pct(80)) %>%
  tab_header(title = md("total number of differentially expressed genes at `adjP<0.05`"))


## Volcano Plots
EnhancedVolcano(DE_Genotype.list$`B6J.LOAD3-B6(Female-4M)`,
                lab = (DE_Genotype.list$`B6J.LOAD3-B6(Female-4M)`$symbol),x = 'log2FoldChange', y = 'padj',legendPosition = 'none',
                title = '4 month Female:LOAD3 vs B6',subtitle = '',
                FCcutoff = 0.0,pCutoff = 0.05,xlim = c(-3, 3))
