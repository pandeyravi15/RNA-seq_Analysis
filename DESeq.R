
# To assess the differentially expressed transcripts. We'll use the `DESeq2` package to perform differential expression analysis of the RNA-seq data.
# custom function for DESeq analysis
DEG.SampleType <- function(rawdata,meta) {
  dseq_res <- data.frame()
  All_res <- data.frame()
  
  dat2 <- as.matrix(rawdata[, rownames(meta)])
  
  ddsHTSeq <- DESeqDataSetFromMatrix(countData = dat2,
                                     colData = meta,
                                     design = ~SampleType)
  
  # Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads in at least 6 separate samples.
  ddsHTSeq <-ddsHTSeq[rowSums(counts(ddsHTSeq) >= 10) >= 6,]
  
  # run the differential expression analysis
  dds <- DESeq(ddsHTSeq, parallel = TRUE)
  res <- results(dds, alpha = 0.05)
  summary(res)
  
  # add gene symbol and ENTREZID to result table
  res$symbol <- map_function.df(res, "ENSEMBL", "SYMBOL")
  res$EntrezGene <- map_function.df(res, "ENSEMBL", "ENTREZID")
  
  # re-arrange resuls table
  All_res <<- as.data.frame(res[, c(7:8, 1:6)])
  
}

# LOAD clean and formatted RNA-seq count data from previous step
load("data/ProcessedData_Brain_Transcriptomics.RData")

# add group column to metadata
metadata <- metadata %>% 
  mutate(Group = paste0(.$Genotype,"-",.$Sex,"-",.$Age,"M")) %>%
  mutate(Genotype = factor(Genotype, levels = c('WT','5XFAD')))

# set up the comparison table to perform DEA
comparisons <-  data.frame(control=c("WT-F-4M" ,  "WT-M-4M" , "WT-F-6M" ,  "WT-M-6M" , "WT-F-12M"  ,  "WT-M-12M"),
                           case=c("5XFAD-F-4M" ,  "5XFAD-M-4M" , "5XFAD-F-6M" ,  "5XFAD-M-6M" , "5XFAD-F-12M"  ,  "5XFAD-M-12M"))

# initiate an empty data frame and list to store results from DEA for each comaprion
DE_Genotype.list <- list()
DE_Genotype.df <- data.frame()

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


## Volcano Plots for differentially expressed genes
for (i in 1:length(DE_Genotype.list))
{
  dat <- DE_Genotype.list[[i]]
  
  print(
    EnhancedVolcano(dat,
                    lab = (dat$symbol),x = 'log2FoldChange', y = 'padj',legendPosition = 'none',
                    title = names(DE_Genotype.list)[i],subtitle = '',
                    FCcutoff = 0.0,pCutoff = 0.05,xlim = c(-5, 5))
  )
}

#save(DE_Genotype.list,DE_Genotype.df,file="~/results/DESeq_Results_Transcripotmics.RData")

