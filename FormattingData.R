#Let's analyze the 5xFAD RNA-seq expression data we explored yesterday. Specifically, we want to know which genes are differentially expressed at each age as a result of the transgenes that constitute the 5xFAD model.

#The 5xFAD mouse model is a 5x transgenic model consisting of mutatnt human transgenes of the amyloid precursor protein (APP) and presenilin 1 (PSEN1) genes. The specific variants are all causal variants for Familial Alzheimer's Disease (FAD) and include three variants in the APP gene - Swedish (K670N, M671L), Florida (I716V), and London (V717I) - and two in the PSEN1 gene - M146L and L286V. The expression of both transgenes is under control of the neural-specific elements of the mouse Thy1 promoter, which drives overexpression of the transgenes in the brain. More information about this generation and maintenance of this strain can be obtained from the [JAX strain catalog](https://www.jax.org/strain/008730).  



### Reading the count data
gene_counts <- read_tsv("../data/rsem.merged.gene_counts.tsv") %>% dplyr::select(-"transcript_id(s)") 
gene_tpm <- read_tsv("../data/rsem.merged.gene_tpm.tsv") %>% dplyr::select(-"transcript_id(s)") 

# Let’s check how many gene_ids are NOT from the mouse genome 
gene_counts[,1:3] %>% 
  filter(!str_detect(gene_id, "MUS"))

#### Accounting for transgenes

#The LOAD3 model has two copies each of the APOE and MAPT genes - one endogenous mouse gene, and the orthologous human transgene. The RNA-seq data was assessed using a custom transcriptome definition that included the sequences of both the mouse and human versions of each gene.  
#Ultimately we are going to sum the counts from both orthologous genes (human APOE and mouse Apoe; human MAPT and mouse Mapt). But first, let's look at the expression of each of these genes in the different groups. To start we'll filter the counts down to just those four relevant gene IDs and join the counts up with the covariates to explore the expression of these genes.  

tg.counts <- gene_counts %>%
  filter(gene_id %in% c("ENSG00000130203","ENSMUSG00000002985",
                        "ENSG00000186868","ENSMUSG00000018411",
                        "ENSG00000264589","ENSG00000279685")) %>% 
  pivot_longer(.,cols = -"gene_id",names_to = "Names",values_to="counts") %>% 
  mutate(Names=as.integer(Names)) %>%
  left_join(meta ,by="Names")

head(tg.counts)

# Let's do a little data housekeeping:

tg.counts <- tg.counts %>% 
  mutate(
    Age = factor(Age, levels = c(4,12))
  )

# add gene symbols
tg.counts <- tg.counts %>% 
  mutate(
    symbol = case_when(
      gene_id == "ENSG00000130203" ~ "Human APOE",
      gene_id == "ENSG00000186868" ~ "Human MAPT",
      gene_id == "ENSMUSG00000002985" ~ "Mouse Apoe",
      gene_id == "ENSMUSG00000018411" ~ "Mouse Mapt",
      gene_id == "ENSG00000264589" ~ "Human MAPT-AS1",
      gene_id == "ENSG00000279685" ~ "Human MAPT-IT1"
    )
  )


#Okay, now let's plot the counts for each gene across the samples:

ggplot(tg.counts, aes(x=Genotype, y=counts, color=Age, shape = Sex)) +
  geom_boxplot() + 
  #geom_point(position=position_jitterdodge())+
  facet_wrap(~symbol, scales = 'free')+
  theme_bw()

#The human transgenes all have a counts of zero in the B6 animals (where the transgenes are absent), while the endogenous mouse genes are expressed relatively consistently across both groups. 
#`Note`: There appears to be one B6 sample, 110507, got sequenced twice.So, either I. skip it.

#Let's combine the expression of corresponding human and mouse genes by summing the expression and saving the summed expression as expression of mouse genes, respectively to match with gene names in control mice. 

# move the gene_id column to rownames, to enable summing across rows
gene_counts <- gene_counts %>% column_to_rownames("gene_id") 

#merge mouse and human APOE gene raw count
gene_counts[rownames(gene_counts) %in% "ENSMUSG00000002985",] <- 
  gene_counts[rownames(gene_counts) %in% "ENSMUSG00000002985",] + 
  gene_counts[rownames(gene_counts) %in% "ENSG00000130203",]

gene_counts <- gene_counts[!rownames(gene_counts) %in% c("ENSG00000130203"),]

#merge mouse and human MAapt gene raw count
gene_counts[rownames(gene_counts) %in% "ENSMUSG00000018411",] <- 
  gene_counts[rownames(gene_counts) %in% "ENSMUSG00000018411",] + 
  gene_counts[rownames(gene_counts) %in% "ENSG00000186868",]

gene_counts <- gene_counts[!rownames(gene_counts) %in% c("ENSG00000186868"),]

# dropping counts from human MAPT-IS1 and MAPT-IT1
gene_counts <- gene_counts[!rownames(gene_counts) %in% c("ENSG00000264589"),]
gene_counts <- gene_counts[!rownames(gene_counts) %in% c("ENSG00000279685"),]

# We can confirm that the human genes are now absent from the counts table:

gene_counts[,1:6] %>% filter(!str_detect(rownames(.), "MUS"))


# convert counts to integer format for DESeq analysis
df_int <- gene_counts %>% mutate(across(everything(), as.integer))
rawcountdata <- df_int[,colnames(df_int) %in% meta$Names]
