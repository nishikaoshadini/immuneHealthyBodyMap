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
library(CellChat)

## from http://tr.im/hH5A

# get_metadata_local(cache_directory = "/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.1.sqlite") |>  filter(.cell == "GCACTAATCCTGGCTT_TSP1_endopancreas_1") |> select(.cell, file_id, lineage_1)


differential_composition_age_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_absolute_FALSE.rds"
data_for_immune_proportion_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/input_absolute.rds"
proportions_age_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_absolute_FALSE_proportion_adjusted.rds"

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

## from http://tr.im/hH5A
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)

data_for_immune_proportion_relative_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/input_relative.rds"
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



  
  
 # ABSOLUTE

# # Save data for third party
# data_for_immune_proportion |>
#
#   # Drop only-immune organs
#   filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone")) |>
#   mutate(is_immune = as.character(is_immune)) |>
#
#   # Mutate days
#   mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>
#   filter(development_stage!="unknown") |>
#   saveRDS("~/PostDoc/sccomp_dev/dev/data_age_absolute_for_third_party.rds")


# Track of immune system in life
differential_composition_age = readRDS(differential_composition_age_absolute_file)

# ggplot() +
#   geom_abline(intercept = -0.281, slope = 0.276)

# generate points
line_age_absolute_mean =
  seq(-3, 3, by = 0.01) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_effect)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_effect)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))

line_age_absolute_lower =
  seq(-3, 3, by = 0.1) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_lower)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_lower)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))

line_age_absolute_upper =
  seq(-3, 3, by = 0.1) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_upper)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_upper)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))



proportions_age_absolute = readRDS(proportions_age_absolute_file)

life_stages = tibble(
  start = c(0, 2, 5, 13, 20, 40, 60),
  end = c(1, 4, 12, 19, 39, 59,100)
) |>
  mutate(stage = c("Infant", "Toddler", "Child", "Teen", "Adult", "Middle age", "Senior"))


# Add life stages
rectangles_age =
  line_age_absolute_mean |>
  mutate(stage = case_when(
  (x_corrected / 365) |> between(0,1) ~ "Infant",
  (x_corrected / 365) |> between(1,4) ~ "Toddler",
  (x_corrected / 365) |> between(4,12) ~ "Child",
  (x_corrected / 365) |> between(12,19) ~ "Teen",
  (x_corrected / 365) |> between(19, 39) ~ "Adult",
  (x_corrected / 365) |> between(39, 59) ~ "Middle age",
  (x_corrected / 365) |> between(59,85) ~ "Senior"
)) |>
  left_join(
    tibble(
      start = c(0, 1, 4, 12, 19, 39, 59),
      end = c(1, 4, 12, 19, 39, 59,85)
    ) |>
      mutate(stage = c("Infant", "Toddler", "Child", "Teen", "Adult", "Middle age", "Senior"))
  ) |>
  with_groups(stage, ~ .x |> mutate(mean_proportion = mean(proportion, na.rm=TRUE))) |>
  distinct(stage, start, end, mean_proportion) |>
  mutate(color =  RColorBrewer::brewer.pal(n = 9, name = "Greys") |> head( n()) )


plot_age_absolute =
  proportions_age_absolute |>
  left_join(data_for_immune_proportion |>
              tidybulk::pivot_sample(.sample)) |>

  filter(development_stage != "unknown") |>

  # Filter
  filter(is_immune == "TRUE") |>
  #filter(tissue_harmonised != "blood") |>
  # Fix samples with multiple assays
  unite(".sample", c(.sample , assay), remove = FALSE) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE) |>

  # Plot
  ggplot() +
  geom_rect(
    aes(xmin = start * 365, xmax = end * 365, ymin = 0, ymax = mean_proportion),
    data = rectangles_age,
    fill = rectangles_age |> select(stage, color) |> deframe(),
    alpha = 0.5
    ) +
  geom_point(
    aes(age_days_original, adjusted_proportion, fill = tissue_harmonised),
    shape = 21,
    stroke = 0,
    size = 1
  ) +
  geom_line(aes(x_corrected, proportion), data = line_age_absolute_mean) +
  geom_line(aes(x_corrected, proportion),
            data = line_age_absolute_lower,
            color = "grey") +
  geom_line(aes(x_corrected, proportion),
            data = line_age_absolute_upper,
            color = "grey") +
  facet_wrap( ~ is_immune, ncol = 9) +
  scale_fill_manual(values = tissue_color) +
  scale_y_continuous( labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none") +
  theme_multipanel

# Color for human heatmap

# Plot age absolute organ cell type
age_absolute_organ_cell_type =
	differential_composition_age |>

  # Find stats of random effect with groups
  test_contrasts(
    contrasts =
      differential_composition_age |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  filter(is_immune == "TRUE")


colors_palette_for_organ_abundance =
	age_absolute_organ_cell_type |>
  select(parameter, c_effect) |> mutate(color = circlize::colorRamp2(
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


# Significance global statistics
count_significance_age_immune_load =
	differential_composition_age |>
	test_contrasts(test_composition_above_logit_fold_change = 0.1) |>
	filter(parameter=="age_days") |>
	filter(is_immune=="TRUE") |>
	count(c_FDR<0.05)


# RELATIVE
differential_composition_age_relative_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_relative_FALSE.rds"
proportions_age_relative_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_relative_FALSE_proportion_adjusted.rds"


differential_composition_age_relative =  readRDS(differential_composition_age_relative_file)

proportions_age_relative = readRDS(proportions_age_relative_file)

# differential_composition_age_relative |>
# 	test_contrasts(test_composition_above_logit_fold_change = 0.2) |>
# 	arrange(desc(abs(c_effect))) |> filter(covariate == "age_days") |> filter(c_FDR<0.05)


line_age_relative_mean =

  differential_composition_age_relative |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  nest(data = -cell_type_harmonised) |>
  mutate(x = list(seq(-3, 3, by = 0.1))) |>
  mutate(y = map2(data, x, ~ {
    .y *
      (.x |> filter(parameter == "age_days") |> pull(c_effect)) +
      (.x |> filter(parameter == "(Intercept)") |> pull(c_effect))

  })) |>
  dplyr::select(-data) |>
  unnest(c(x, y)) |>
  with_groups(x, ~ .x |> mutate(proportion = softmax(y))) |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))


# Volcano relative
library(scales)
library(ggplot2)

S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

volcano_relative = 
 differential_composition_age_relative |>
  
  filter(parameter == "age_days") |> 
  mutate(naive_experienced = case_when(
    cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
    cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
  )) |> 
  mutate(significant = c_FDR<0.05) |> 
  mutate(cell_type_harmonised = case_when(c_FDR<0.05~cell_type_harmonised)) |> 
  ggplot(aes(c_effect, c_FDR)) + 
  geom_point(aes(color=naive_experienced, size=significant)) +
  ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
  scale_x_continuous(trans="S_sqrt") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_size_discrete(range = c(0, 0.5)) +
  theme_multipanel

# Plot age relative cell types
plot_age_relative =
  proportions_age_relative |>
  
  inner_join(
    differential_composition_age_relative |>

      filter(parameter == "age_days") |>
      filter(c_FDR<0.05) |>
      distinct(cell_type_harmonised, c_effect)
  ) |>
  
  left_join(data_for_immune_proportion_relative |>
              tidybulk::pivot_sample(.sample)) |>

  filter(development_stage != "unknown") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  #filter(tissue_harmonised != "blood") |>
  # Fix samples with multiple assays
  unite(".sample", c(.sample , assay), remove = FALSE) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE)  |>

  # Filter
  filter(cell_type_harmonised != "immune_unclassified")  |> 

  # Relevel
  arrange(c_effect) %>%
  mutate(cell_type_harmonised = factor(cell_type_harmonised, levels = unique(.$cell_type_harmonised))) |> 
  
  ggplot(aes(age_days_original, adjusted_proportion)) +
  geom_point(
    aes(fill = tissue_harmonised),
    shape = 21,
    stroke = 0,
    size = 0.4
  ) +
  geom_line(
    aes(x_corrected, proportion, color = significant),
    data = line_age_relative_mean |>

      # Join statistics
      inner_join(
        differential_composition_age_relative |>
          filter(parameter == "age_days") |>
          mutate(significant = c_FDR < 0.05) |>
          filter(cell_type_harmonised != "immune_unclassified") |>
          dplyr::select(cell_type_harmonised, significant) |> 
          filter(significant)
      )
  ) +
  facet_wrap( ~ cell_type_harmonised, ncol = 3) +
  scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  scale_fill_manual(values = tissue_color) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none", color = "none") +
  theme_multipanel


# Plot age relative organ cell type
age_relative_organ_cell_type =

  differential_composition_age_relative |>
  test_contrasts(
    contrasts =
      differential_composition_age_relative |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + {parameter}") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( )
  )



plot_age_relative_by_organ =

  age_relative_organ_cell_type |>

  filter(c_FDR<0.01) |>
  arrange(desc(abs(c_effect))) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  add_count(cell_type_harmonised) |>
  arrange(parameter, desc(n)) |>
  ggplot(aes( fct_reorder(cell_type_harmonised, desc(n)), parameter)) +
  geom_tile(aes(fill = c_effect)) +
  scale_fill_distiller(palette="Spectral") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=20, hjust = 1, vjust = 1))




df_heatmap_age_relative_organ_cell_type =

  differential_composition_age_relative |>

  # Find stats of random effect with groups
  test_contrasts(
    contrasts =
      differential_composition_age_relative |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( ),
    test_composition_above_logit_fold_change = 0.4
  )  |>

  filter(cell_type_harmonised != "immune_unclassified") |>
  add_count(cell_type_harmonised) |>
  arrange(parameter, desc(n)) |>

  rename(tissue = parameter) |>
  rename(cell_type = cell_type_harmonised) |>

  # To be fixed in the model
  mutate(is_treg = cell_type =="treg") |> 
  nest(data = -is_treg) |> 
  mutate(data = map2(
    data, is_treg,
    ~ {
      if(.y) .x |> mutate(c_effect = c_effect/7 )
      else(.x)
    }
  )) |> 
  unnest(data) |> 
  
  # Cell type abundance
  with_groups(cell_type, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(cell_type_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # Filter for visualisation
  filter(!cell_type %in% c("non_immune", "immune_unclassified")) |>

  # Tissue diversity
  with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # First rank
  with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

  # # Cap
  # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>
  mutate(Difference = c_effect) |>
  
  rename(`Mean diff` = cell_type_mean_change) |>
  mutate(`Mean diff tissue` = -tissue_mean_change) |>
  mutate(cell_type = cell_type |> str_replace("macrophage", "macro")) |>
  mutate(tissue = tissue |> str_replace_all("_", " ")) |>

  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>

  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>

  # Shorten names
  mutate(cell_type = cell_type |> 
           str_replace("megakaryocytes", "mega") |> 
           str_remove("phage") |> 
           str_replace("th1/th17", "th1/17") |> 
           str_replace("mono", "mn") |> 
           str_replace("tcm", "cm") |> 
           str_replace("cd4 th", "Th") |> 
           str_replace("memory", "mem") |> 
           str_replace("naive", "nv") |> 
           str_replace("terminal effector cd4 t", "cd4 eff") |> 
           str_remove("cyte") 
  ) |> 
  mutate(tissue = tissue |> 
           str_replace("intestine", "int") |> 
           str_replace("large", "lrg") |> 
           str_replace("small", "sml") |> 
           str_replace("node", "nd") |> 
           str_replace("prostate", "prost")
  ) |> 
  
  # Order
  mutate(tissue = fct_reorder(tissue, `Mean diff tissue`)) |>
  mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))


plot_heatmap_age_relative_organ_cell_type =

  df_heatmap_age_relative_organ_cell_type |>
  
  # Heatmap
  heatmap(
    tissue, cell_type, Difference,
    # palette_value = circlize::colorRamp2(
    #   seq(-3, 3, length.out = 11),
    #   RColorBrewer::brewer.pal(11, "Spectral")
    # ),
    palette_value = circlize::colorRamp2(
      seq(3, -3, length.out = 11),
      RColorBrewer::brewer.pal(11, "RdBu")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    show_heatmap_legend = FALSE
  ) |>

  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  annotation_tile(
    tissue, show_legend = FALSE,
    palette =
      df_heatmap_age_relative_organ_cell_type |>
      distinct(tissue, tissue_color) |>
      arrange(tissue) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    cell_type, show_legend = FALSE,
    palette =
      df_heatmap_age_relative_organ_cell_type |>
      distinct(cell_type, cell_type_color)  |>
      arrange(cell_type) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    count_tissue, show_legend = FALSE,
    size = unit(0.2, "cm"),
    palette = c( "white", "black")
  ) |>
  layer_point((c_lower * c_upper)>0)


plot_heatmap_age_relative_organ_cell_type |>
  save_pdf(
    filename = "~/PostDoc/CuratedAtlasQueryR/dev/plot_heatmap_age_relative_organ_cell_type.pdf",
    width = 80*1.5, height = 60*1.5, units = "mm"
  )

# Color body
df_heatmap_age_relative_organ_cell_type |>
  distinct(tissue, `Mean diff tissue`) |> 
  arrange(`Mean diff tissue`) |> 
  mutate(
    color =  circlize::colorRamp2(c(0, 5, 10, 12), viridis::viridis(4))(-`Mean diff tissue`)
  ) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>  pull(color) |>
  scales::show_col()



count_significance_age_cell_type =
  differential_composition_age_relative |>
  test_contrasts(test_composition_above_logit_fold_change = 0.4) |>
  filter(parameter=="age_days") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)

count_significance_age_cell_type_tissue =
  df_heatmap_age_relative_organ_cell_type |>
  count(c_FDR<0.05)



rm(differential_composition_age_relative , differential_composition_age )
gc()

# plot_significance_overall =
#   count_significance_age_immune_load |>
#   mutate(name = c(
#     "count_significance_age_immune_load"
#   )) |>
#   bind_rows(
#     count_significance_age_immune_load_tissue |>
#       mutate(name = c(
#         "count_significance_age_immune_load_tissue"
#       ))
#     ) |>
#   bind_rows(count_significance_age_cell_type |>
#               mutate(name = c(
#
#                 "count_significance_age_cell_type"
#               ))) |>
#   bind_rows(count_significance_age_cell_type_tissue |>
#               mutate(name = c(
#
#                 "count_significance_age_cell_type_tissue"
#               ))) |>
#
#   bind_rows(count_significance_sex_immune_load |>
#               mutate(name = c(
#
#                 "count_significance_sex_immune_load"
#               ))) |>
#   bind_rows(count_significance_sex_immune_load_tissue |>
#               mutate(name = c(
#
#                 "count_significance_sex_immune_load_tissue"
#               ))) |>
#   bind_rows(count_significance_sex_cell_type |>
#               mutate(name = c(
#
#                 "count_significance_sex_cell_type"
#               ))) |>
#   bind_rows(count_significance_sex_cell_type_tissue |>
#               mutate(name = c(
#
#                 "count_significance_sex_cell_type_tissue"
#               ))) |>
#
#   bind_rows(count_significance_ethnicity_immune_load |>
#               mutate(name = c(
#
#                 "count_significance_ethnicity_immune_load"
#               ))) |>
#   bind_rows(count_significance_ethnicity_immune_load_tissue |>
#               mutate(name = c(
#
#                 "count_significance_ethnicity_immune_load_tissue"
#               ))) |>
#   bind_rows(count_significance_ethnicity_cell_type |>
#               mutate(name = c(
#
#                 "count_significance_ethnicity_cell_type"
#               ))) |>
#   bind_rows(count_significance_ethnicity_cell_type_tissue |>
#               mutate(name = c(
#
#                 "count_significance_ethnicity_cell_type_tissue"
#               ))) |>
#     tidyr::extract(name, c("factor", "variable", "resolution"), "count_significance_([a-zA-Z]+)_([a-zA-Z]+_[a-zA-Z]+)_?(.*)", remove = FALSE) |>
#     mutate(resolution = if_else(resolution == "", "overall", resolution)) |>
#     unite("xlab", c(variable, resolution), remove = FALSE) |>
#     mutate(xlab = xlab |> fct_relevel(c("immune_load_overall", "immune_load_tissue", "cell_type_overall", "cell_type_tissue"))) |>
#     with_groups(name, ~ .x |> mutate(sum_n = sum(n))) |>
#     mutate(proportion = n/sum_n) |>
#     mutate(factor = factor |> str_to_sentence()) |>
#     ggplot(aes(xlab, proportion, fill=`c_FDR < 0.05`)) +
#     geom_bar(stat = "identity")+
#     geom_text(aes(y = 0.5, label = sum_n), size = 2.5, angle=90) +
#     facet_wrap( ~ factor,  nrow=1) +
#     scale_fill_manual(values = c("FALSE"="grey", "TRUE"="#D5C711")) +
#     ylab("Proportion of significant tests") +
#     xlab("Hypotheses") +
#     theme_multipanel +
#     theme(axis.text.x = element_text(angle=20, hjust = 1, vjust = 1))



# job::job({ readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_relative_4.rds") |> remove_unwanted_variation(~ age_days) })

# Residency
plot_trends_residency =
  readRDS("~/PostDoc/CuratedAtlasQueryR/dev/residency_data.rds") |>
  
  filter(cell_type_harmonised |> str_detect("naive", negate = T)) |>
  separate(cell_type_harmonised, "cell_type_harmonised", sep = " ") |>
  
  # Filter significant
  inner_join(
    readRDS("~/PostDoc/CuratedAtlasQueryR/dev/residency_estimates_random_effects.rds") |> 
      filter((`l-90% CI` * `u-90% CI`) > 0) |> 
      filter(Marker == "residency_axel") |> 
      
      # Filter only with a lot of data
      filter(tissue_harmonised %in% (
        df_heatmap_age_relative_organ_cell_type |> 
          distinct(tissue) |> 
          pull(tissue) |>
          as.character()
      )) |> 
      distinct(cell_type_harmonised, tissue_harmonised, Marker) 
  ) |> 
  
  # Slope
  nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, Marker)) |>
  mutate(slope = map_dbl(data, ~ lm(mean_TotalScore ~ age_days,data = .x)$coeff[2] )) |>
  unnest(data) |> 
  
  # Filter only tissue to plot
  filter(tissue_harmonised %in% c("adipose", "bone", "intestine small", "spleen", "heart", "uterus")) |> 
  mutate(tissue_harmonised = tissue_harmonised |> fct_relevel(c("adipose", "bone", "intestine small", "spleen", "heart", "uterus"))) |> 
  
  ggplot(aes(age_days, mean_TotalScore)) +
  geom_point(size=0.2) +
  geom_smooth(aes(group=glue("{file_id} {cell_type_harmonised}"), color=slope>0), method="lm", se=F, size=0.3) +
  #facet_grid(tissue_harmonised ~ facet) +
  facet_wrap(~ tissue_harmonised, nrow=1, ) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC")) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  xlab("Years") +
  ylab("Mean residency score") +
  ylim(-1, NA) + 
  theme_multipanel




plot_first_line =
  (
    plot_spacer() | #plot_significance_overall |
      ((plot_age_absolute / plot_spacer()) + plot_layout(heights  = c(3,2))) |
      plot_spacer()
  ) +
  plot_layout(widths = c(46,96, 36),  guides = 'collect') 


plot_second_line =
  (
    volcano_relative | 
      plot_age_relative  |
      plot_spacer() |
      wrap_heatmap(plot_heatmap_age_relative_organ_cell_type, padding = unit(c(-40, 0, -10, -30), "points" ))
  ) +
  plot_layout( widths = c(0.1, 0.2, 0.23, 0.37), guides = 'collect') 

plot_third_line = plot_trends_residency | plot_spacer()


# Plotting
p =
  (
    plot_first_line /
      plot_second_line /
      plot_third_line 
  ) +
  plot_layout(guides = 'collect', heights = c(57, 57, 30, 57))  &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )


ggsave(
  "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/figure_age.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 200 ,
  limitsize = FALSE
)



