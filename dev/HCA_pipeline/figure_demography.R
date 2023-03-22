library(tidyverse)
library(forcats)
library(CuratedAtlasQueryR)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
library(tidyHeatmap)
library(ComplexHeatmap)
library(scales)
library(ggplot2)

## from http://tr.im/hH5A

# get_metadata_local(cache_directory = "/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.1.sqlite") |>  filter(.cell == "GCACTAATCCTGGCTT_TSP1_endopancreas_1") |> select(.cell, file_id, lineage_1)
result_directory = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1"

softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }

  exp(x - logsumexp(x))
}

dropLeadingZero <-
	function(l) {
		stringr::str_replace(l, '0(?=.)', '')
	}
S_sqrt <- function(x) {
	sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
	x ^ 2 * sign(x)
}
S_sqrt_trans <-
	function()
		scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

clean_names = function(x){
  x |>  mutate(
    tissue_harmonised =
      tissue_harmonised |>
      str_remove("tissue_harmonised") |>
      str_replace_all("_", " ") |>
      str_replace("gland", "gld") |>
      str_replace("node", "nd") |>
      str_replace("skeletal", "sk")
  )
}

circle_plot = function(res) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }
  softmax <- function (x) {
    exp(x - logsumexp(x))
  }
  
  res_relative_for_plot =
    res |>
    filter(!parameter |> str_detect("group___")) |>
    
    # Cell type abundance
    with_groups(tissue_harmonised, ~ .x |>  mutate(proportion = softmax(c_effect))) |>
    with_groups(cell_type_harmonised,
                ~ .x |>  mutate(cell_type_mean_abundance = mean(proportion))) |>
    
    # Filter for visualisation
    filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>
    
    # Tissue diversity
    with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>
    
    # First rank
    with_groups(cell_type_harmonised,
                ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>
    
    # Cap
    mutate(c_effect = c_effect |> pmax(-5) |> pmin(5))
  
  inter_type_diversity_plot =
    res_relative_for_plot |>
    distinct(inter_type_diversity, tissue_harmonised) |>
    ggplot(aes(
      inter_type_diversity,
      fct_reorder(tissue_harmonised, inter_type_diversity)
    )) +
    geom_bar(stat = "identity") +
    scale_x_reverse() +
    xlab("Diversity") +
    ylab("Tissue") +
    theme_multipanel
  
  cell_type_mean_abundance_plot =
    res_relative_for_plot |>
    distinct(cell_type_mean_abundance, cell_type_harmonised) |>
    ggplot(aes(
      fct_reorder(
        cell_type_harmonised,
        dplyr::desc(cell_type_mean_abundance)
      ),
      cell_type_mean_abundance
    )) +
    geom_bar(stat = "identity") +
    scale_y_continuous(position = "right") +
    theme_multipanel +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
    )
  
  circle_plot =
    res_relative_for_plot |>
    arrange(rank == 1) |>
    ggplot() +
    geom_point(aes(
      fct_reorder(
        cell_type_harmonised,
        dplyr::desc(cell_type_mean_abundance)
      ),
      fct_reorder(tissue_harmonised, inter_type_diversity) ,
      fill = rank,
      size = c_effect,
      stroke = rank == 1
    ),
    shape = 21) +
    scale_size_continuous(range = c(0.5, 3)) +
    xlab("Cell type") +
    theme_multipanel +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  
  plot_spacer() +
    cell_type_mean_abundance_plot +
    inter_type_diversity_plot +
    circle_plot +
    plot_layout(guides = 'collect',
                height = c(1, 4),
                width = c(1, 8)) &
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          legend.position = "bottom")
  
}


## from http://tr.im/hH5A
data_for_immune_proportion_absolute_file = glue("{result_directory}/input_absolute.rds")
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)

data_for_immune_proportion_relative_file = glue("{result_directory}/input_relative.rds")
data_for_immune_proportion_relative = readRDS(data_for_immune_proportion_relative_file)

tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()


source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")

cell_type_color = 
  data_for_immune_proportion |> 
  pull(cell_type_harmonised) |> 
  unique() |> 
  get_cell_type_color()

names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")


# SEX
differential_composition_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE.rds")
differential_composition_sex_relative_file = glue("{result_directory}/sex_relative_FALSE.rds")
proportions_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE_proportion_adjusted.rds")

differential_composition_sex_relative = 
  readRDS(differential_composition_sex_relative_file)  

differential_composition_sex_absolute = 
  readRDS(differential_composition_sex_absolute_file)


proportions_sex_absolute_adjusted = readRDS(proportions_sex_absolute_file)

draws_abundance = 
  readRDS(differential_composition_sex_absolute_file)   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE")

draws_variability = 
  readRDS(differential_composition_sex_absolute_file)   |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE") |> filter(parameter=="sexmale")

plot_sex_absolute_1D =
  tibble(
    Variability = draws_variability |> pull(.value),
    Abundance = draws_abundance |> pull(.value)
  ) |>
  tidybulk::as_matrix() |>
  bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
  coord_flip() +
  xlab("Effect male immune cellularity") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=90, hjust = 0.5))


# plot_sex_absolute =
#   proportions_sex_absolute_adjusted |>
#   left_join(data_for_immune_proportion |> distinct(.sample, sex, age_days, tissue_harmonised)) |>
#   rename(proportion = adjusted_proportion) |>
#   bind_rows(
#     differential_composition_sex_absolute  |> 
#       sccomp_replicate(~ sex, ~ sex) |>
#       left_join(data_for_immune_proportion |> distinct(.sample, sex, age_days)) |>
#       mutate(tissue_harmonised = "Underlying dist.") |>
#       rename(proportion = generated_proportions)
#   ) |>
#   filter(is_immune=="TRUE") |>
#   
#   ggplot(aes(y = tissue_harmonised, x = proportion)) +
#   ggridges::geom_density_ridges(aes(fill = tissue_harmonised)) +
#   facet_grid(sex ~.) +
#   scale_fill_manual(values = tissue_color) +
#   theme_multipanel


sex_absolute_organ_tissue =
  differential_composition_sex_absolute |> 
  test_contrasts(
    contrasts =
      differential_composition_sex_absolute |>
      filter(parameter |> str_detect("___sex")) |>
      distinct(parameter) |>
      mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      filter(contrast |> str_detect("_female", negate = TRUE)) |> 
      deframe( ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  filter(is_immune == "TRUE")

colors_palette_for_organ_abundance =
  sex_absolute_organ_tissue |>
  select(parameter, c_effect) |> 
  mutate(color = circlize::colorRamp2(
    seq(1.45,-1.45, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )(c_effect)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>
  pull(color) |>
  scales::show_col(	cex_label = 0.5	)



plot_sex_absolute_organ_boxoplot_adjusted =
  proportions_sex_absolute_adjusted |>
  left_join(
    data_for_immune_proportion |>
      distinct(.sample, tissue_harmonised, ethnicity, sex,tissue, file_id)
  ) |>
  inner_join(
    sex_absolute_organ_tissue |> 
      filter(c_FDR<0.07) |> 
      separate(parameter, c("tissue_harmonised", "sex"), sep="_") |>  
      distinct(tissue_harmonised)
  ) |> 
  
  mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
  filter(is_immune =="TRUE") |>
  ggplot(aes(sex, adjusted_proportion )) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
  geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
  facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
  guides(color = "none") +
  ylab("Adjusted proportion") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))





# Sex relative

# Volcano relative
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

volcano_relative_sex = 
  differential_composition_sex_relative |>
  
  test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
  filter(cell_type_harmonised != "non_immune") |> 
  
  filter(parameter == "sexmale") |> 
  mutate(naive_experienced = case_when(
    cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
    cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
  )) |> 
  mutate(significant = c_FDR<0.05) |> 
  mutate(cell_type_harmonised = case_when(c_FDR<0.05~cell_type_harmonised)) |> 
  ggplot(aes(c_effect, c_FDR)) + 
  
  geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
  
  geom_point(aes(color=naive_experienced, size=significant)) +
  ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
  scale_x_continuous(trans="S_sqrt") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_size_discrete(range = c(0, 0.5)) +
  guides(size="none", color="none") +
  theme_multipanel

# df_heatmap_sex_relative_organ_cell_type =
#   
#   differential_composition_sex_relative |> 
# 
#   # Find stats of random effect with groups
#   test_contrasts(
#     contrasts =
#       differential_composition_sex_relative |>
#       filter(parameter |> str_detect("_male___sex")) |>
#       distinct(parameter) |>
#       mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
#       deframe( ),
#     test_composition_above_logit_fold_change = 0.1
#   )  |>
#   
#   filter(cell_type_harmonised != "immune_unclassified") |>
#   add_count(cell_type_harmonised) |>
#   arrange(parameter, desc(n)) |>
#   
#   rename(tissue = parameter) |>
#   rename(cell_type = cell_type_harmonised) |>
#   
#   # To be fixed in the model
#   mutate(is_treg = cell_type =="treg") |> 
#   nest(data = -is_treg) |> 
#   mutate(data = map2(
#     data, is_treg,
#     ~ {
#       if(.y) .x |> mutate(c_effect = c_effect/7 )
#       else(.x)
#     }
#   )) |> 
#   unnest(data) |> 
#   
#   # Cell type abundance
#   with_groups(cell_type, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(cell_type_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>
#   
#   # Filter for visualisation
#   filter(!cell_type %in% c("non_immune", "immune_unclassified")) |>
#   
#   # Tissue diversity
#   with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>
#   
#   # First rank
#   with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>
#   
#   # # Cap
#   # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>
#   mutate(Difference = c_effect) |>
#   
#   rename(`Mean diff` = cell_type_mean_change) |>
#   mutate(`Mean diff tissue` = -tissue_mean_change) |>
#   mutate(cell_type = cell_type |> str_replace("macrophage", "macro")) |>
#   mutate(tissue = tissue |> str_replace_all("_", " ")) |>
#   
#   # Color
#   left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
#   left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>
#   
#   # Counts
#   left_join(
#     data_for_immune_proportion_relative |>
#       count(tissue_harmonised, name = "count_tissue") |>
#       rename(tissue = tissue_harmonised) |>
#       mutate(count_tissue = log(count_tissue))
#   ) |>
#   
#   # Shorten names
#   mutate(cell_type = cell_type |> 
#            str_replace("megakaryocytes", "mega") |> 
#            str_remove("phage") |> 
#            str_replace("th1/th17", "th1/17") |> 
#            str_replace("mono", "mn") |> 
#            str_replace("tcm", "cm") |> 
#            str_replace("cd4 th", "Th") |> 
#            str_replace("memory", "mem") |> 
#            str_replace("naive", "nv") |> 
#            str_replace("terminal effector cd4 t", "cd4 eff") |> 
#            str_remove("cyte") 
#   ) |> 
#   mutate(tissue = tissue |> 
#            str_replace("intestine", "int") |> 
#            str_replace("large", "lrg") |> 
#            str_replace("small", "sml") |> 
#            str_replace("node", "nd") |> 
#            str_replace("prostate", "prost")
#   ) |> 
#   
#   # Order
#   mutate(tissue = fct_reorder(tissue, `Mean diff tissue`)) |>
#   mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))
# 
# plot_heatmap_sex_relative_organ_cell_type =
#   
#   df_heatmap_sex_relative_organ_cell_type |>
#   
#   # Heatmap
#   heatmap(
#     tissue, cell_type, Difference,
# 
#     palette_value = circlize::colorRamp2(
#       seq(1, -1, length.out = 11),
#       RColorBrewer::brewer.pal(11, "RdBu")
#     ),
#     cluster_rows = FALSE,
#     cluster_columns = FALSE,
#     row_names_gp = gpar(fontsize = 6),
#     column_names_gp = gpar(fontsize = 6),
#     column_title_gp = gpar(fontsize = 0),
#     row_title_gp = gpar(fontsize = 0),
#     show_heatmap_legend = FALSE
#   ) |>
#   
#   annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
#   annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
#   annotation_tile(
#     tissue, show_legend = FALSE,
#     palette =
#       df_heatmap_sex_relative_organ_cell_type |>
#       distinct(tissue, tissue_color) |>
#       arrange(tissue) |>
#       deframe(),
#     size = unit(0.2, "cm")
#   ) |>
#   annotation_tile(
#     cell_type, show_legend = FALSE,
#     palette =
#       df_heatmap_sex_relative_organ_cell_type |>
#       distinct(cell_type, cell_type_color)  |>
#       arrange(cell_type) |>
#       deframe(),
#     size = unit(0.2, "cm")
#   ) |>
#   annotation_tile(
#     count_tissue, show_legend = FALSE,
#     size = unit(0.2, "cm"),
#     palette = c( "white", "black")
#   ) |>
#   layer_point((c_lower * c_upper)>0)


# differential_composition_sex_relative |>
#   remove_unwanted_variation(~ sex, ~ sex) |>
#   saveRDS("~/PostDoc/HCAquery/dev/proportions_sex_relative_adjusted.rds")


# proportions_sex_absolute_adjusted |>
#   filter(cell_type_harmonised == "cd14 mono") |>
#   left_join(data_for_immune_proportion_relative |> distinct(.sample, sex, age_days, tissue_harmonised)) |>
# 
#   rename(proportion = adjusted_proportion) |>
# 
#   ggplot(aes(y = proportion, x = sex)) +
#   geom_violin() +
#   scale_fill_manual(values = tissue_color) +
#   theme_multipanel


# Ethnicicy

differential_composition_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE.rds")
proportions_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE_proportion_adjusted.rds")

differential_composition_ethnicity_relative_file = glue("{result_directory}/ethnicity_relative_FALSE.rds")

differential_composition_ethnicity_relative = readRDS(differential_composition_ethnicity_relative_file)
#differential_composition_ethnicity_relative = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1_BACKUP/ethnicity_relative_FALSE.rds")

differential_composition_ethnicity_absolute = readRDS(differential_composition_ethnicity_absolute_file)

proportions_ethnicity_tissue_absolute_adjusted = readRDS(proportions_ethnicity_absolute_file)

  
data_for_ethinicity_absolute_plot =
  differential_composition_ethnicity_absolute |>
  test_contrasts(
    c(
      African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
      Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
      European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
      Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
    )
  ) |>
  
  filter(parameter %in% c("African",  "Hispanic", "European", "Chinese")) |>
  # rowwise() |>
  # mutate(
  #   proportion_mean = softmax(c(c_effect,-c_effect))[1],
  #   proportion_lower = softmax(c(c_lower,-c_lower))[1],
  #   proportion_upper = softmax(c(c_upper,-c_upper))[1]
  # ) |>
  # ungroup() |>
  filter(is_immune == "TRUE")


plot_diff_abundance =
  data_for_ethinicity_absolute_plot |>
  ggplot(aes(c_effect, fct_reorder(parameter, desc(c_effect)))) +
  geom_vline(xintercept = 0, color="grey", linetype = "dashed") +
  geom_linerange(aes(xmin = c_lower, xmax = c_upper), lwd = 0.2,) +
  geom_point(size = 0.5, color = "red") +
  xlab("Effect compared to mean of other ethnicictes") +
  ylab("Ethnicity") +
  theme_multipanel


plot_diff_variability =
  data_for_ethinicity_absolute_plot |>
  ggplot(aes(v_effect, fct_reorder(parameter, desc(c_effect)))) +
  geom_vline(xintercept = 0, color="grey", linetype = "dashed") +
  geom_linerange(aes(xmin = v_lower, xmax = v_upper), lwd = 0.2,) +
  geom_point(size = 0.5, color = "red") +
  xlab("Variability log-fold change") +
  ylab("Ethnicity") +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )


draws_abundance = 
  differential_composition_ethnicity_absolute   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(
    African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
    Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
    European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
    Asian = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
  )) |> 
  filter(is_immune=="TRUE")

draws_variability = 
  differential_composition_ethnicity_absolute  |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(
    African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
    Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
    European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
    Asian = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
  )) |> 
  filter(is_immune=="TRUE") |> 
  filter(parameter %in% c("Hispanic", "African", "European", "Asian")) 

plot_ethnicity_absolute_1D =
  left_join(
    draws_variability |> select(.value, parameter, .draw) |> rename(Variability = .value),
    draws_abundance |> select(.value, parameter, .draw) |> rename(Abundance = .value),
  ) |>
  select(-.draw) |> 
  select(parameter, everything()) |> 
  nest(data = -parameter) |> 
  mutate(plot = map2(
    data, parameter,
    ~ .x |> 
      tidybulk::as_matrix() |>
      bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
      coord_flip() +
      xlab("Effect male immune cellularity") +
      xlim(-0.75, 0.75) +
      theme_multipanel +
      theme(axis.text.x = element_text(angle=90, hjust = 0.5)) +
      ggtitle(.y)
  )) |> 
  rownames_to_column() |> 
  mutate(plot = map2(
    plot, rowname, 
    ~ if(.y>1) .x + 
      theme(
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()
      ) else .x)) |> 
  pull(plot) |> 
  wrap_plots(nrow=1)
  



ethnicity_absolute_organ =
  differential_composition_ethnicity_absolute |>
  
  # Find stats of random effect with groups
  test_contrasts(
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
      select(tissue_harmonised, contrast) |> 
      separate(tissue_harmonised, c("tissue_harmonised", "ethnicity"), sep="_") |> 
      add_count(tissue_harmonised) |> mutate(contrast = glue("({contrast})")) |>  
      pivot_wider(names_from = ethnicity, values_from = contrast) |> 
      
      # Filter where there are at least 3 ethnicities
      mutate(n_minus_1 = n-1) |> 
      filter(n_minus_1>1) |>
      
      # Build contrasts
      rowwise() |> 
      mutate(
        other_Chinese = c(European, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        other_European = c(Chinese, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        other_African = c(European, Chinese, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        `other_Hispanic or Latin American` = c(European, African, Chinese) |> str_subset(".+") |> str_c(collapse=" + "),
      ) |> 
      mutate(
        contrast_Chinese = glue("1/{n_minus_1}*({other_Chinese}) - {Chinese}"),
        contrast_European = glue("1/{n_minus_1}*({other_European}) - {European}"),
        contrast_African = glue("1/{n_minus_1}*({other_African}) - {African}"),
        contrast_Hispanic= glue("1/{n_minus_1}*({`other_Hispanic or Latin American`}) - {`Hispanic or Latin American`}"),
      ) |> 
      rename(
        this_Chinese = Chinese,
        this_European = European,
        this_African = African,
        this_Hispanic = `Hispanic or Latin American`
      ) |> 
      select(tissue_harmonised, starts_with(c("contrast_", "this_"))) |>
      pivot_longer(cols = starts_with(c("contrast_", "this_")), names_to = c("type", "ethnicity"), names_sep = "_") |> 
      pivot_wider(names_from = type, values_from = value) |>
      filter(!is.na(this)) |> 
      select(tissue_harmonised, ethnicity, contrast) |> 
      unite("tissue_harmonised", c(tissue_harmonised, ethnicity)) |> 
      deframe( )
  )

# # CI Absolute
# ethnicity_absolute_organ |>
#   filter(is_immune == "TRUE") |>
#   tidyr::extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>
#   
#   nest(data = -tissue) |> 
#   filter(map_dbl(data, ~ .x |> pull(c_FDR) |> min()) <0.05) |> 
#   unnest(data) |> 
#   # # Select significant
#   # inner_join(
#   #   differential_composition_ethnicity_absolute_one_vs_all |>
#   #     extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>
#   #     distinct(tissue)
#   # ) |>
#   # extract(parameter, "tissue", "(.+)_", remove = FALSE) |>
#   ggplot(aes(c_effect, fct_reorder(parameter, dplyr::desc(c_effect)))) +
#   geom_linerange(aes(xmin = c_lower, xmax = c_upper), lwd = 0.2,) +
#   geom_point(size = 0.5, color = "red") +
#   facet_wrap(~tissue, scale="free_y") +
#   xlab("Immune proportion") +
#   ylab("Ethnicity") +
#   theme_multipanel

# Body map
ethnicity_absolute_organ |>
  filter(is_immune == "TRUE") |>
  select(parameter, c_effect, c_FDR) |> 
  separate(parameter, c("tissue", "ethnicity"), remove = FALSE, sep="_") |> 
  mutate(color = circlize::colorRamp2(
    seq(2,-2, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )(c_effect)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>
  
  # group
  nest(data = -ethnicity) |> 
  mutate(dummy = map(data, ~ scales::show_col(.x$color,	cex_label = 0.5	))) |> 
  pull(data)


plot_ethnicity_absolute_organ_boxoplot_adjusted =
  proportions_ethnicity_tissue_absolute_adjusted |>
  left_join(
    data_for_immune_proportion |>
      distinct(.sample, tissue_harmonised, ethnicity, tissue, file_id)
  ) |>
  inner_join(
    ethnicity_absolute_organ |> 
      filter(c_FDR<0.05) |>  
      separate(parameter, c("tissue", "ethnicity"), remove = FALSE, sep="_") |> 
      distinct(tissue) |> rename(tissue_harmonised = tissue)
  ) |>
  mutate(ethnicity = ethnicity |> str_remove(" or Latin American") ) |>
  mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
  filter(is_immune =="TRUE") |>
  ggplot(aes(ethnicity, adjusted_proportion )) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
  geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
  facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
  guides(color = "none") +
  ylab("Adjusted proportion") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))




# Relative
data_for_ethinicity_relative_plot =
  differential_composition_ethnicity_relative |>
  test_contrasts(
    c(
      African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
      Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
      European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
      Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
    ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  
  filter(parameter %in% c("African",  "Hispanic", "European", "Chinese"))

plot_ethinicy_bubble =
  data_for_ethinicity_relative_plot |>
  mutate(tissue_harmonised = parameter) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  circle_plot() +
  scale_fill_viridis_c(direction = -1)


# significant_cell_types =
#   differential_composition_ethnicity_relative |>
#   test_contrasts(
#     c(
#       African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
#       Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
#       European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
#       Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
#     )
#   ) |>
#   
#   filter(parameter %in% c("African",  "Hispanic", "European", "Chinese")) |>
#   nest(data = -cell_type_harmonised) |>
#   mutate(sd_effect = map_dbl(data, ~ sd(.x$c_effect))) |>
#   filter(cell_type_harmonised != "immune_unclassified") |>
#   arrange(desc(sd_effect)) |>
#   head(6) |>
#   pull(cell_type_harmonised)


S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

volcano_relative_ethnicity = 
  data_for_ethinicity_relative_plot |>
  filter(cell_type_harmonised != "non_immune") |> 
  mutate(naive_experienced = case_when(
    cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
    cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
  )) |>
  mutate(significant = c_FDR<0.05) |> 
  mutate(cell_type_harmonised = case_when(c_FDR<1e-2~cell_type_harmonised)) |> 
  ggplot(aes(c_effect, c_FDR)) + 
  geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
  
  geom_point(aes(color=naive_experienced, size=significant)) +
  ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
  facet_wrap(~parameter, nrow=1) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
  scale_x_continuous(trans="S_sqrt") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_size_discrete(range = c(0, 0.5)) +
  theme_multipanel



rm(differential_composition_ethnicity_relative , differential_composition_ethnicity_absolute )
gc()


plot = 
  ((
    ((
      plot_sex_absolute_1D |
        plot_spacer() |
        plot_sex_absolute_organ_boxoplot_adjusted |
        volcano_relative_sex
    ) + plot_layout(width = c(12, 22, 12, 30))) /
      
      plot_spacer() /
      
      (plot_ethnicity_absolute_1D  + plot_layout(width = c(2,2,2,2)) ) /
      
      volcano_relative_ethnicity 
      
  ) + plot_layout( height = c(55, 18, 20, 50))
) |
  ((
    plot_spacer() / plot_ethnicity_absolute_organ_boxoplot_adjusted / plot_spacer()
  )  + plot_layout( height = c(85, 50, 38)) ) +
  plot_layout(width = c(84, 97)) &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )




ggsave(
  glue("~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/figure_demography.pdf"),
  plot = plot,
  units = c("mm"),
  width = 183 ,
  height = 180 ,
  limitsize = FALSE
)

# job::job({ readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_4.rds") |> remove_unwanted_variation(~ age_days) })

# Plotting
p = 
 (
   ( plot_sex | plot_spacer() )  +
    plot_layout(guides = 'collect', widths =  c(45, 55))
  )  &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )






differential_composition_sex_relative |> 
  filter(c_FDR<0.05) |> 
  filter(parameter != "sexunknown") |> 
  filter(factor=="sex") |> 
  print(n=99)


