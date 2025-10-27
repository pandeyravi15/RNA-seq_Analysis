### Overview of Human Consensus RNA-Seq Coexpression Modules
#[Wan, et al.](https://doi.org/10.1016/j.celrep.2020.107908) performed multi
#method co-expression network analysis  followed by differential analysis and
#found 30 co-expression modules related LOAD pathology from human cohort study.
#Among the 30 aggregate co-expression modules, five consensus clusters have been
#described by Wan, et al. These consensus clusters consist of a subset of
#modules which are associated with similar AD related changes across the
#multiple studies and brain regions.These AMP-AD co-expression modules are very
#useful for making comparisons between animal models and the human cohorts.

#In order to use these modules to make the comparisons, we'll need to download
#data pertaining to the 30 co-expression modules. These data are available from
#the Synapse data repository
#([syn11932957](https://www.synapse.org/Synapse:syn11932957)); let's download
#and take a closer look at the data.

query <- syn$tableQuery("SELECT * FROM syn11932957")
module_table <- read_csv(query$filepath)

#We can visualize this as bar plot using ggplot2 package. 
ggplot(module_table,aes(y=Module)) + 
  geom_bar() + 
  theme_bw() 

### Mouse orthologs of Human module genes  
#In the module table we've downloaded we have human ENSEMBL ids and human gene
#symbols. In order to compare between human and mouse models, we need to
#identify the corresponding (i.e. orthologous) mouse gene IDs. We are going to
#add the gene IDs of orthologous genes in mouse to the corresponding human genes
#in the module table. Orthology mapping can be tricky, but  Wan *et al* have
#already identified mouse orthologs for each of the human genes using the HGNC
#Comparison of Orthology Predictions
#([HCOP](https://www.genenames.org/tools/hcop/)) tool.

#We are going to read that table from Synapse ([syn17010253](https://doi.org/10.7303/syn17010253.1)).

mouse.human.ortho <- syn$get("syn17010253")$path %>% read_tsv()
#mouse.human.ortho <- read_tsv("human_mouse_hcop_fifteen_column_20160523.txt")

module_table$Mouse_gene_symbol <-
  mouse.human.ortho$mouse_symbol[
    match(module_table$GeneID,mouse.human.ortho$human_ensembl_gene)
  ]

#Some genes don't have identified orthologs. Also there's some redundant
#information. Let's only keep the columns of interest and rows that contain a
#mouse ortholog mapping:

ampad_modules <- module_table %>%
  distinct(tissue = brainRegion, module = Module, gene = GeneID, Mouse_gene_symbol) %>%
  filter(Mouse_gene_symbol != "")


### Reading differential expression result of human data from meta-analysis

#Now we know the genes that are in each AMP-AD co-expression cluster, along with
#the ID of the corresponding orthologous genes in mouse. We'll also need to know
#how the expression of these genes change in AD.

#We'll download the results from differential expression analysis of reprocessed
#AMP-AD RNA-seq data (all 7 brain regions). Log fold change values for human
#transcripts can be obtained from Synapse
#([syn14237651](https://www.synapse.org/#!Synapse:syn14237651)).
                                                                                                                                                                                                                      
ampad_modules_raw <- read_tsv(syn$get("syn14237651")$path)

# The AMP-AD data has been processed many ways and using different models and comparisons. 
# For our analyses we'll use data from the "Diagnosis" model and comparisons
#between AD cases vs controls. We'll filter the table for these conditions and
#only keep the three columns we'll need: `Tissue`, `Gene` and `logFC`.

ampad_fc <- ampad_modules_raw %>%
  as_tibble() %>%
  filter(Model == "Diagnosis", Comparison == "AD-CONTROL") %>%
  dplyr::select(tissue = Tissue, gene = ensembl_gene_id, ampad_fc = logFC)

#Next, we will combine the fold change table (`ampad_fc`) and
#module table (`ampad_modules`).

ampad_modules_fc <- ampad_modules %>%
  inner_join(ampad_fc, by = c("gene", "tissue")) %>% 
  dplyr::select(symbol = Mouse_gene_symbol, module, ampad_fc) 

# We will use the data in `ampad_modules_fc` to compare with log fold change data measured in mouse models.

#### Preparing module information for correlation plots 

mod <- c("TCXblue","PHGyellow","IFGyellow","DLPFCblue","CBEturquoise","STGblue","PHGturquoise","IFGturquoise","TCXturquoise","FPturquoise","IFGbrown","STGbrown","DLPFCyellow","TCXgreen","FPyellow","CBEyellow","PHGbrown","DLPFCbrown","STGyellow","PHGgreen","CBEbrown","TCXyellow","IFGblue","FPblue","FPbrown","CBEblue","DLPFCturquoise","TCXbrown","STGturquoise","PHGblue")

cluster_a <- tibble(
  module = c("TCXblue", "PHGyellow", "IFGyellow"),
  cluster = "Consensus Cluster A (ECM organization)",
  cluster_label = "Consensus Cluster A\n(ECM organization)"
)

cluster_b <- tibble(
  module = c("DLPFCblue", "CBEturquoise", "STGblue", "PHGturquoise", "IFGturquoise", "TCXturquoise", "FPturquoise"),
  cluster = "Consensus Cluster B (Immune system)",
  cluster_label = "Consensus Cluster B\n(Immune system)"
)

cluster_c <- tibble(
  module = c("IFGbrown", "STGbrown", "DLPFCyellow", "TCXgreen", "FPyellow", "CBEyellow", "PHGbrown"),
  cluster = "Consensus Cluster C (Neuronal system)",
  cluster_label = "Consensus Cluster C\n(Neuronal system)"
)

cluster_d <- tibble(
  module = c("DLPFCbrown", "STGyellow", "PHGgreen", "CBEbrown", "TCXyellow", "IFGblue", "FPblue"),
  cluster = "Consensus Cluster D (Cell Cycle, NMD)",
  cluster_label = "Consensus Cluster D\n(Cell Cycle, NMD)"
)

cluster_e <- tibble(
  module = c("FPbrown", "CBEblue", "DLPFCturquoise", "TCXbrown", "STGturquoise", "PHGblue"),
  cluster = "Consensus Cluster E (Organelle Biogensis, Cellular stress response)",
  cluster_label = "Consensus Cluster E\n(Organelle Biogenesis,\nCellular stress response)"
)

module_clusters <- cluster_a %>%
  bind_rows(cluster_b) %>%
  bind_rows(cluster_c) %>%
  bind_rows(cluster_d) %>%
  bind_rows(cluster_e) %>%
  mutate(cluster_label = fct_inorder(cluster_label))

mod <- module_clusters$module

save(ampad_modules_fc,module_clusters,mod, file= here::here("results","AMPAD_Module_Data.RData"))

#There are two approaches that we commonly use to compute correlation between
#mouse data and human AD data. Both approaches allow us to assess directional
#coherence between gene expression for genes in AMP-AD modules and the effects
#of genetic manipulations in mice. For this session we are going to use the
#first approach; we'll return to approach #2 later in the week.

# [1] Compare change in expression in Human AD cases vs controls with change in expression in mouse models for each gene in a given module:  
#  + LogFC(h) = log fold change in expression of human AD patients compared to control patients. 
#  + LogFC(m) = log fold change in expression of mouse AD models compared to control mouse models.

cor.test(LogFC(h), LogFC(m))
  
#Read the results saved after differential expression analysis
load("data/DESeq_Results_Transcripotmics.RData")

# we'll combine both mouse and human `ampad_modules_fc` log fold change data
# sets for all genes and compute correlation coefficients using the `cor.test` function
mouse_vs_ampad_fc <- DE_Genotype.df %>%
  inner_join(ampad_modules_fc,
             by = c('symbol'),
             multiple = "all") %>%
  dplyr::select(module, group, symbol, log2FoldChange, ampad_fc) %>%
  group_by(module, group) %>%
  nest(data = c(symbol, log2FoldChange, ampad_fc)) %>%
  mutate(
    cor_test = map(data, ~ cor.test(.x[["log2FoldChange"]], .x[["ampad_fc"]], method = "pearson")),
    estimate = map_dbl(cor_test, "estimate"),
    p_value = map_dbl(cor_test, "p.value"),
    coherent = map_dbl(data, ~ sum((.x[["log2FoldChange"]]*.x[["ampad_fc"]])>0)/nrow(.x)),
    anticoherent = map_dbl(data, ~ sum((.x[["log2FoldChange"]]*.x[["ampad_fc"]])<0)/nrow(.x))
  ) %>%
  ungroup() %>%
  dplyr::select(-cor_test)

#### Annotate correlation table to prepare for visualization
#We'll add a column that flags whether the correlation is significant or
#not, and we'll add in the information about which consensus cluster each module
#belongs to:
  
mouse_vs_ampad <- mouse_vs_ampad_fc %>%
  mutate(significant = p_value < 0.05) %>%
  left_join(module_clusters, by = "module") %>%
  dplyr::select(cluster,cluster_label,module,group,
                correlation = estimate,p_value,significant)

#### Create a data frame to use as input for plotting the results

data_for_plot <- mouse_vs_ampad %>%
  arrange(cluster) %>%
  mutate(
    module = factor(module, levels = mod),
    group = factor(group,levels=c("5XFAD-B6(F-4M)","5XFAD-B6(M-4M)","5XFAD-B6(F-6M)","5XFAD-B6(M-6M)", "5XFAD-B6(F-12M)","5XFAD-B6(M-12M)"))
  ) 

### Visualizing the Correlation plot
data <- data_for_plot
range(data_for_plot$correlation)

p <- ggplot() +
  
  # the AMP-AD modules will be along the x-axis, the mouse models will be along the y-axis
  geom_tile(data = data, aes(x = .data$module, y = .data$group), colour = "black", fill = "white") +
  
  # each tile of the grid will be filled with a circle where the fill and size correspond to the correlation coefficient
  geom_point(data = data, aes(x = .data$module, y = .data$group, 
                              colour = .data$correlation, size = abs(.data$correlation))) +
  
  # we'll draw a box arround significant correlations
  geom_point(data = dplyr::filter(data, .data$significant), 
             aes(x = .data$module, y = .data$group, colour = .data$correlation),
             color="black",shape=0,size=9) +
  
  # plot the x-axis on the top of the plot, set the parameters of the scales
  scale_x_discrete(position = "top") +
  scale_size(guide = "none", limits = c(0, 0.4)) +
  scale_color_gradient2(limits = c(-0.5, 0.5), low = "#85070C", high = "#164B6E", 
                        name = "Correlation", guide = guide_colorbar(ticks = FALSE)) +
  
  # remove axis labels
  labs(x = NULL, y = NULL) +
  
  # facet the plot based on age range for the mice (rows) and consensus cluster (columns)
  facet_grid( cols = dplyr::vars(.data$cluster_label), 
              scales = "free", space = "free",switch = 'y') +
  
  # specify how different aspects of the plot will look
  theme(
    strip.text.x = element_text(size = 10,colour = c("black")),
    strip.text.y.left = element_text(angle = 0,size = 12),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, size = 12),
    axis.text.y = element_text(angle = 0, size = 12),
    panel.background = element_blank(),
    plot.title = element_text(angle = 0, vjust = -54, hjust = 0.03,size=12,face="bold"),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    legend.position = "right"
  )

fills <- c("darkorange3","chartreuse3","deepskyblue2","turquoise","deeppink2")
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
ggdraw(g)


#In the above plot, categories along the x-axis are the 30 AMP-AD co-expression
#modules grouped into 5 consensus clusters, while the categories along the
#y-axis show the different groupings of mouse models tested (split by age).
#Positive correlations are shown in blue and negative correlations in red. 
#Color intensity and size of the circles are proportional to the correlation
#coefficient.  Black squares around dots represent significant correlation at
#p-value=0.05 and non-significant correlations are left blank.

#########################################################################################################
### 2. Compare Human AD expression changes to mouse genetic effects for each gene in a given module:  ###
#########################################################################################################

#  + h = human gene expression (Log2 RNA-seq Fold Change AD/control)
#  + β = mouse gene expression effect from linear regression model (Log2 RNA-seq TPM)

cor.test(LogFC(h), β)

# Here, we are using log normalized TPM count data to assess the effect size of each factors
mydat <- as.data.frame(mydat_TPM) %>% 
  mutate(Gene=mapIds(org.Mm.eg.db,keys = rownames(.),column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first"),.before =1) %>%
  pivot_longer(cols=-Gene,names_to="Names",values_to="value") %>% 
  na.omit(.) %>%
  mutate(Names = as.integer(Names))

mydat_with_metadata <- mydat %>% left_join(metadata,by="Names") %>% 
  dplyr::select(Gene,Names,Sex,Age,Genotype,Age,value) %>% 
  mutate(Sex = factor(Sex,levels=c("F","M")))  %>% 
  mutate(FAD = ifelse(Genotype %in% c("5XFAD"),1,0))  %>%
  mutate(Age=factor(Age, levels = c(4,6,12)))

# run the multiple linear regression model
lms.case <- mydat_with_metadata %>% split(.$Gene) %>% map(~lm(value ~ Sex + Age + FAD, data=.x)) %>% map(summary) 

# extract the effect size and corresponding p-values
effects <- as.data.frame(bind_rows( lapply(lms.case, function(x) x$coefficients[,"Estimate"] ))) %>% 
  mutate(names = names(lms.case)) %>% column_to_rownames(.,var="names") %>% 
  tibble::rownames_to_column(.,"Gene") %>% dplyr::select(-"(Intercept)") %>% 
  pivot_longer(cols=-Gene,names_to="Variant",values_to="effect_size") %>%
  mutate(Variant = case_when(
    Variant == "SexM" ~ "Sex (Male)",
    Variant == "Age6" ~ "Age_6M",
    Variant == "Age12" ~ "Age_12M",
    Variant == "FAD" ~ "5XFAD"
  ))

pvals <- as.data.frame(bind_rows( lapply(lms.case, function(x) x$coefficients[,"Pr(>|t|)"] ))) %>% 
  mutate(names = names(lms.case)) %>% column_to_rownames(.,var="names") %>% 
  tibble::rownames_to_column(.,"Gene") %>% dplyr::select(-"(Intercept)") %>% 
  pivot_longer(cols=-Gene,names_to="Variant",values_to="pval") %>%
  mutate(Variant = case_when(
    Variant == "SexM" ~ "Sex (Male)",
    Variant == "Age6" ~ "Age_6M",
    Variant == "Age12" ~ "Age_12M",
    Variant == "FAD" ~ "5XFAD"
  ))

# join effect size and p-value for each factor for each gene in a single daat frame.
effects_pvals <- effects %>% left_join(pvals,by=c("Gene","Variant")) %>% mutate(FDR = p.adjust(pval,method = "fdr"))

# Created a function to perform correlation analysis and generate data frame for plotting like previous section.
corr_function_lm <- function(data, human_data)
{
  mouse_vs_ampad <- data  %>%
    inner_join(human_data, by = c("Gene")) %>%
    group_by(module, Variant) %>%
    nest(data = c(Gene, effect_size, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["effect_size"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  # Flag for significant results, add cluster information to modules
  mouse_vs_ampad <- mouse_vs_ampad %>%
    mutate(significant = p_value < 0.05) %>%
    left_join(module_clusters, by = "module") %>%
    dplyr::select(
      cluster, cluster_label, module,Variant,
      correlation = estimate, p_value,significant
    )
  
  # Create a version of the data for plotting - clean up naming, order factors, etc
  data_for_plot.all <- mouse_vs_ampad %>%
    arrange(cluster) %>%
    mutate(
      Variant = factor(Variant, levels = ordered.variant),
      Variant = fct_rev(Variant),
      module = factor(module, levels = mod),
    )
}

ordered.variant <- c("Sex (Male)","Age_6M","Age_12M","5XFAD")
data_for_plot <- corr_function_lm(effects,ampad_modules_fc %>% rename("Gene"="symbol")) %>% 
  mutate(Background = ifelse(Variant %in% c("Sex (Male)"),"Female",ifelse(Variant %in% c("5XFAD"),"B6","Age_4M"))) %>%
  mutate(Background=factor(Background,levels = c("Female","Age_4M","B6")))
range(data_for_plot$correlation)


### Visualizing the Correlation plot
range(data_for_plot$correlation)

# Function to create correlation plot
variant_corrplot_colored <- function(data, ran) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data, 
      ggplot2::aes(x = .data$module, y = .data$Variant),
      colour = "black", fill = "white"
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
        colour = .data$correlation),
      color = "black",
      shape = 0,
      size = 9) +
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
        face = "bold"),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) + theme(legend.title = element_text(size=13),legend.text = element_text(size=12))
  
  fills <- c("darkorange3","chartreuse3","deepskyblue2","turquoise","deeppink2")
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-t', g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  ggdraw(g)
  
}


variant_corrplot_colored(data_for_plot,0.4)
