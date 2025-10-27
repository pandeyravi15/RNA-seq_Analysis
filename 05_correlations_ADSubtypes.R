#Compute correlation with AD subtypes in ROSMAP, Mayo and MSBB cohort
#identified by [Nikhil et al.](https://pubmed.ncbi.nlm.nih.gov/32492070/) and
#five AD subtypes in MSBB cohort identified by [Neff et
#al.](https://pubmed.ncbi.nlm.nih.gov/33523961/). Nikhil's subtype's were
#annotated as inflammatory(ROSMAP_subtypeA, Mayo_subtypeA & B, MSBB_subtypeA)
#and non-inflammatory subtypes. Neff's subtypes were classified into three
#larger classes: typical AD (subtype C1 & C2), intermediate (subtype B1 & B2),
#or atypical AD (subtype A).

#Read the results saved after differential expression analysis
load("data/DESeq_Results_Transcripotmics.RData")

### Correlation with Nikhil's subtypes using logFC values

order.model <- c("5XFAD-B6(F-4M)","5XFAD-B6(M-4M)","5XFAD-B6(F-6M)","5XFAD-B6(M-6M)", "5XFAD-B6(F-12M)","5XFAD-B6(M-12M)")
plot_subtype1 <- subtypecorr_function(DE_Genotype.df) 

names(load2_fc)[names(load2_fc) == 'symbol'] <- 'Gene'
ns_vs_ampad_fc <- DE_Genotype.df %>%
  inner_join(ampad_subtype_fc, by = c("Gene")) %>%
  dplyr::select(subtype, group, Gene, log2FoldChange, ampad_fc) %>%
  group_by(subtype, group) %>%
  nest(data = c(Gene, log2FoldChange, ampad_fc)) %>%
  mutate(
    cor_test = map(data, ~ cor.test(.x[["log2FoldChange"]], .x[["ampad_fc"]], method = "pearson")),
    estimate = map_dbl(cor_test, "estimate"),
    p_value = map_dbl(cor_test, "p.value")
  ) %>%
  ungroup() %>%
  dplyr::select(-cor_test)

nanostring <- ns_vs_ampad_fc %>%
  mutate(significant = p_value < 0.05) %>%
  left_join(module_cohort, by = "subtype") %>%
  dplyr::select(cluster,
                cluster_label,
                subtype,
                group,
                correlation = estimate,
                p_value,
                significant)

nanostring_for_plot_subtype <- nanostring %>%
  arrange(cluster) %>%
  mutate(
    subtype = factor(subtype, levels = subtypes),
    model_sex = group,
  ) %>%
  arrange(model_sex) %>%
  mutate(
    model_sex = factor(model_sex, levels = order.model),
    model_sex = fct_rev(model_sex),
  )



range(plot_subtype1$correlation)
subtype_FC_corrplot(plot_subtype1,0.15)

### Correlation with Neff's subtypes using logFC values
order.model <- c("B6J.LOAD3-B6(Female-4M)" ,  "B6J.LOAD3-B6(Male-4M)" , "B6J.LOAD3-B6(Female-12M)","B6J.LOAD3-B6(Male-12M)")
plot_subtype2 <- Neff_subtypecorr_function(DE_Genotype.df)  
range(plot_subtype2$correlation)
subtype_FC_corrplot(plot_subtype2,0.15)

subtypecorr_function <- function(load2_fc)
{
  names(load2_fc)[names(load2_fc) == 'symbol'] <- 'Gene'
  ns_vs_ampad_fc <- load2_fc %>%
    inner_join(ampad_subtype_fc, by = c("Gene")) %>%
    dplyr::select(subtype, group, Gene, log2FoldChange, ampad_fc) %>%
    group_by(subtype, group) %>%
    nest(data = c(Gene, log2FoldChange, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["log2FoldChange"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < 0.05) %>%
    left_join(module_cohort, by = "subtype") %>%
    dplyr::select(cluster,
                  cluster_label,
                  subtype,
                  group,
                  correlation = estimate,
                  p_value,
                  significant)
  
  nanostring_for_plot_subtype <- nanostring %>%
    arrange(cluster) %>%
    mutate(
      subtype = factor(subtype, levels = subtypes),
      model_sex = group,
    ) %>%
    arrange(model_sex) %>%
    mutate(
      model_sex = factor(model_sex, levels = order.model),
      model_sex = fct_rev(model_sex),
    )
}

Neff_subtypecorr_function <- function(load2_fc)
{
  names(load2_fc)[names(load2_fc) == 'symbol'] <- 'Gene'
  ns_vs_ampad_fc <- load2_fc %>%
    inner_join(neff_subtype_fc , by = c("Gene")) %>%
    dplyr::select(subtype, group, Gene, log2FoldChange, ampad_fc) %>%
    group_by(subtype, group) %>%
    nest(data = c(Gene, log2FoldChange, ampad_fc)) %>%
    mutate(
      cor_test = map(data, ~ cor.test(.x[["log2FoldChange"]], .x[["ampad_fc"]], method = "pearson")),
      estimate = map_dbl(cor_test, "estimate"),
      p_value = map_dbl(cor_test, "p.value")
    ) %>%
    ungroup() %>%
    dplyr::select(-cor_test)
  
  nanostring <- ns_vs_ampad_fc %>%
    mutate(significant = p_value < 0.05) %>%
    #left_join(module_cohort, by = "subtype") %>%
    dplyr::select(subtype, group, correlation = estimate, p_value, significant)
  nanostring$cluster_label <- "MSBB_PHG_Neff"
  
  nanostring_for_plot_subtype <- nanostring %>%
    mutate(#subtype = factor(subtype,levels=subtypes),
      model_sex = group) %>%
    arrange(model_sex) %>%
    mutate(
      model_sex = factor(model_sex, levels = order.model),
      model_sex = fct_rev(model_sex),
    )
}

subtype_FC_corrplot <- function(data, ran) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = data,
      ggplot2::aes(x = .data$subtype, y = .data$model_sex),
      colour = "black",
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data),
      ggplot2::aes(
        x = .data$subtype,
        y = .data$model_sex,
        colour = .data$correlation,
        size = abs(.data$correlation)
      )
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, .data$significant),
      ggplot2::aes(
        x = .data$subtype,
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
    ggplot2::facet_grid(
      #rows = dplyr::vars(.data$group),
      cols = dplyr::vars(.data$cluster_label),
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 11, colour = c("blue")),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 11),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = 11,
        angle = 90,
        hjust = 0
      ),
      axis.text.y = ggplot2::element_text(size = 11, angle = 0),
      panel.background = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}