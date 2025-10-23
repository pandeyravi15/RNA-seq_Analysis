###Let's analyze the 5xFAD RNA-seq expression data we explored yesterday. Specifically, we want to know which genes are differentially expressed at each age as a result of the transgenes that constitute the 5xFAD model.The 5xFAD mouse model is a 5x transgenic model consisting of mutatnt human transgenes of the amyloid precursor protein (APP) and presenilin 1 (PSEN1) genes. The specific variants are all causal variants for Familial Alzheimer's Disease (FAD) and include three variants in the APP gene - Swedish (K670N, M671L), Florida (I716V), and London (V717I) - and two in the PSEN1 gene - M146L and L286V. The expression of both transgenes is under control of the neural-specific elements of the mouse Thy1 promoter, which drives overexpression of the transgenes in the brain. More information about this generation and maintenance of this strain can be obtained from the [JAX strain catalog](https://www.jax.org/strain/008730).  

# Reading metadata
meta <- read.csv("data/metadata_5XFAD_RNASeq_JAX.csv")


# Reading the count data
gene_counts <- read_tsv("data/rnaseq_rsem.merged.gene_counts.tsv") %>% dplyr::select(-"transcript_id(s)") 
gene_tpm <- read_tsv("data/rnaseq_rsem.merged.gene_tpm.tsv") %>% dplyr::select(-"transcript_id(s)") 

# Let’s check how many gene_ids are NOT from the mouse genome 
gene_counts[,1:3] %>% 
  filter(!str_detect(gene_id, "MUS"))

#### Accounting for transgenes
tg.counts <- gene_counts %>%
  filter(gene_id %in% c("ENSG00000080815","ENSMUSG00000019969",
                        "ENSG00000142192","ENSMUSG00000022892")) %>% 
  pivot_longer(.,cols = -"gene_id",names_to = "Names",values_to="counts") %>% 
  mutate(Names=as.integer(Names)) %>%
  left_join(meta ,by=c("Names"))

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
      gene_id == "ENSG00000142192" ~ "Human APP",
      gene_id == "ENSG00000080815" ~ "Human PSEN1",
      gene_id == "ENSMUSG00000022892" ~ "Mouse App",
      gene_id == "ENSMUSG00000019969" ~ "Mouse Psen1"
    )
  )


#Okay, now let's plot the counts for each gene across the samples:

ggplot(tg.counts, aes(x=Genotype, y=counts, color=Age, shape = Sex)) +
  geom_boxplot() + 
  #geom_point(position=position_jitterdodge())+
  facet_wrap(~symbol, scales = 'free')+
  theme_bw()

#The human transgenes all have a counts of zero in the WT animals (where the transgenes are absent), while the endogenous mouse genes are expressed relatively consistently across both groups. 

#Let's combine the expression of corresponding human and mouse genes by summing the expression and saving the summed expression as expression of mouse genes, respectively to match with gene names in control mice.
# move the gene_id column to rownames, to enable summing across rows
counts <- gene_counts %>% column_to_rownames("gene_id") 

#merge mouse and human APOE gene raw count
counts[rownames(counts) %in% "ENSMUSG00000022892",] <-  counts[rownames(counts) %in% "ENSMUSG00000022892",] + counts[rownames(counts) %in% "ENSG00000142192",]

#merge mouse and human PS1 gene raw count
counts[rownames(counts) %in% "ENSMUSG00000019969",] <- counts[rownames(counts) %in% "ENSMUSG00000019969",] + counts[rownames(counts) %in% "ENSG00000080815",]

counts <- counts[!rownames(counts) %in% c("ENSG00000080815","ENSG00000142192"),]

# We can confirm that the human genes are now absent from the counts table:
counts[,1:6] %>% filter(!str_detect(rownames(.), "MUS"))


# convert counts to integer format for DESeq analysis
df_int <- counts %>% mutate(across(everything(), as.integer))
rawcountdata <- df_int[,colnames(df_int) %in% meta$Names]


###### Let's account for transgenes for gene tpm count data as well 
tg.counts <- gene_tpm %>%
  filter(gene_id %in% c("ENSG00000080815","ENSMUSG00000019969",
                        "ENSG00000142192","ENSMUSG00000022892")) %>% 
  pivot_longer(.,cols = -"gene_id",names_to = "Names",values_to="tpm") %>% 
  mutate(Names=as.integer(Names)) %>%
  left_join(meta ,by="Names")

tg.counts <- tg.counts %>% 
  mutate(
    Age = factor(Age, levels = c(4,12))
  )

# add gene symbols
tg.counts <- tg.counts %>% 
  mutate(
    symbol = case_when(
      gene_id == "ENSG00000142192" ~ "Human APP",
      gene_id == "ENSG00000080815" ~ "Human PSEN1",
      gene_id == "ENSMUSG00000022892" ~ "Mouse App",
      gene_id == "ENSMUSG00000019969" ~ "Mouse Psen1"
    )
  )


#Okay, now let's plot the counts for each gene across the samples:
ggplot(tg.counts, aes(x=Genotype, y=tpm, color=Sex, shape = Age)) +
  geom_boxplot() + 
  geom_point(position=position_jitterdodge())+
  facet_wrap(~symbol, scales = 'free')+
  theme_bw()


#The human transgenes all have a counts of zero in the B6 animals (where the transgenes are absent), while the endogenous mouse genes are expressed relatively consistently across both groups. 

#Let's combine the expression of corresponding human and mouse genes by summing the expression and saving the summed expression as expression of mouse genes, respectively to match with gene names in control mice. 

# move the gene_id column to rownames, to enable summing across rows
tpm <- gene_tpm %>% column_to_rownames("gene_id") 

#merge mouse and human APOE gene raw count
tpm[rownames(tpm) %in% "ENSMUSG00000022892",] <- 
  tpm[rownames(tpm) %in% "ENSMUSG00000022892",] + 
  tpm[rownames(tpm) %in% "ENSG00000142192",]

#merge mouse and human PS1 gene raw count
tpm[rownames(tpm) %in% "ENSMUSG00000019969",] <- 
  tpm[rownames(tpm) %in% "ENSMUSG00000019969",] + 
  tpm[rownames(tpm) %in% "ENSG00000080815",]

# dropping counts from human APP and PSEN1
tpm <- tpm[!rownames(tpm) %in% c("ENSG00000080815","ENSG00000142192"),]


# We can confirm that the human genes are now absent from the counts table:
tpm[,1:6] %>% filter(!str_detect(rownames(.), "MUS"))

#Here we perform a minimal pre-filtering to keep only rows that have at least 0.1  in gene_tpm count data in at least 6 separate samples.
TPMcount.filter <- tpm[rowSums((tpm > 0.1))>5, ]
rawdata1 <- TPMcount.filter[,colnames(TPMcount.filter) %in% meta$Names]
mydat_TPM <- as.data.frame(log2(rawdata1 + 1))



# We will also validate sex of samples. Let's plot the counts for female specific `Xist` and male-specific `Ddx3y` genes across the samples:
 tg.counts <- tpm %>% rownames_to_column(.,var="gene_id") %>%
  filter(gene_id %in% c("ENSMUSG00000086503","ENSMUSG00000069045")) %>% 
  pivot_longer(.,cols = -"gene_id",names_to = "Names",values_to="tpm") %>% 
  mutate(Names=as.integer(Names)) %>%
  left_join(meta ,by="Names")


tg.counts <- tg.counts %>% 
  mutate(
    Age = factor(Age, levels = c(4,12))
  )

tg.counts <- tg.counts %>% 
  mutate(
    symbol = case_when(
      gene_id == "ENSMUSG00000086503" ~ "Mouse Xist",
      gene_id == "ENSMUSG00000069045" ~ "Mouse Ddx3y",
    )
  )

ggplot(tg.counts, aes(x=Sex, y=tpm, color=Age, shape = Genotype)) +
  geom_boxplot() + 
  geom_point(position=position_jitterdodge())+
  facet_wrap(~symbol, scales = 'free')+
  theme_bw()


# Now, we have formatted the count and metadata for downstream analyses. We will save the data
metadata <- as.data.frame(meta)
rownames(metadata) <- metadata$Names

#save(metadata,rawcountdata,mydat_TPM,file="data/ProcessedData_Brain_Transcriptomics.RData")



