library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
library(ggforce)
library(ggpubr)
library(tidyHeatmap)
library(ComplexHeatmap)
library(tidybulk)

# Read arguments
args = commandArgs(trailingOnly = TRUE)
relative_assay =  args[[1]]
cell_metadata_with_harmonised_annotation = args[[2]]
sql_lite_DB = args[[3]]
output_rds = args[[4]]
output_pdf = args[[5]]

## from http://tr.im/hH5A

# Calculate softmax from an array of reals
softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }

  exp(x - logsumexp(x))
}

# This function shorten names of cell types for visualisation purposes
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


# Get the dataset for PCA with uncertainty
data_for_plot_1 =
  readRDS(cell_metadata_with_harmonised_annotation) |>

  left_join(
    get_metadata(sql_lite_DB) |>
      dplyr::select(.cell, is_primary_data.y, name, cell_type, file_id, assay) |>
      as_tibble()
  ) |>

  # Count samples
  distinct(.sample, tissue_harmonised, file_id, assay) |>
  add_count(tissue_harmonised, name = "Sample count") |>

  # Add colours
  nest(data = -assay) |>
  arrange(assay) |>
  mutate( color = RColorBrewer::brewer.pal(7,"Blues") |> c(dittoSeq::dittoColors()[-2][1:7]) ) |>
  unnest(data)


#------------------------------#
# Analyses of immune cellularity proportion of immune cells in a tissue
#------------------------------#

# Read results
relative_assay_results = readRDS(relative_assay)

# Set contrasts that will be used several times
contrasts = c(
  "assay10x_3_v3 - assay10x_3_v2" ,
  "assaysci_RNA_seq - assay10x_3_v2" ,
  "assaymicrowell_seq - assay10x_3_v2" ,
  "assay10x_5_v2 - assay10x_3_v2" ,
  "assay10x_5_v1 - assay10x_3_v2" ,
  "assaySmart_seq2  - assay10x_3_v2"
)

# Plot Heatmap of the difference in immune composition across technologies
plot_heatmap_relative_assay =
  relative_assay_results |>
  filter(parameter %in% contrasts) |>
  mutate(parameter = parameter |> str_remove_all("assay")) |>
  mutate(parameter = parameter |> str_remove("- 10x_3_v2")) |>
  mutate(tissue_harmonised = parameter) |>

  # Calculate stats
  filter(!parameter |> str_detect("group___")) |>

  # Cell type abundance
  with_groups(cell_type_harmonised, ~ .x |>  mutate(cell_type_mean_change = mean(abs(c_effect)))) |>

  # Filter for visualisation
  filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>

  # Tissue diversity
  with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>

  # First rank
  with_groups(cell_type_harmonised, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

  mutate(Difference = c_effect) |>
  mutate(inter_type_diversity = -inter_type_diversity) |>
  rename(`Mean diff` = cell_type_mean_change) |>
  rename(Diversity =inter_type_diversity ) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

  # Order
  mutate(parameter = fct_reorder(parameter, Diversity)) |>
  mutate(cell_type_harmonised = fct_reorder(cell_type_harmonised, -`Mean diff`)) |>


  # Heatmap
  heatmap(
    parameter, cell_type_harmonised, Difference,
    palette_value = circlize::colorRamp2(
      seq(8, -8, length.out = 11),
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
  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
  annotation_bar(Diversity, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
  layer_point((c_lower * c_upper)>0)


## Gat dataset for PCA plot
# data_for_assay_pca =
# 	relative_assay_results |>
#   sccomp:::get_abundance_contrast_draws(
#     contrasts = c(
#       "assay10x_3_v2",
#       "assay10x_3_v3" ,
#       "assaysci_RNA_seq" ,
#       "assaymicrowell_seq" ,
#       "assay10x_5_v2" ,
#       "assay10x_5_v1" ,
#       "assaySmart_seq2"
#     )
#   )|>
#   unite("sample", c(parameter, .draw), remove = FALSE) |>
#   reduce_dimensions(
#     sample, cell_type_harmonised, .value,
#     method="PCA", action="get", scale=FALSE,
#     transform = identity
#   )
#
# data_for_assay_pca |> saveRDS("~/PostDoc/HCAquery/dev/data_for_assay_pca.rds")

data_for_assay_pca = readRDS("~/PostDoc/HCAquery/dev/data_for_assay_pca.rds")

# Plot PCA with uncertainty 
plot_assay_PCA =
  data_for_assay_pca |>
  mutate(parameter = parameter |> str_remove_all("assay")) |>
  ggplot(aes(PC1, PC2, label = parameter, color = parameter)) +
  stat_density_2d(geom = "polygon", aes(fill = parameter),  alpha = 0.2, bins = 2) +
  geom_point(aes(color = parameter), data = data_for_assay_pca |> with_groups(parameter, ~ .x |> summarise(across(c(PC1, PC2), median) ))) +
  scale_fill_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all(" |-", "_") |> str_remove("'")) |>  deframe() ) +
  scale_color_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all(" |-", "_") |> str_remove("'")) |>  deframe() ) +
  guides(fill="none", color = "none")  +
  theme_multipanel


# Get analysis of differential variability across assays
# This plots how consistent each assay is across cell types
# excluding other sources of variability, including
# Sex, ethnicity, age, tissue and random effects including datasets
res_relative_for_variability_plot =
	relative_assay_results |>
  test_contrasts(
    contrasts = c(
      "assay10x_3_v2",
      "assay10x_3_v3" ,
      "assaysci_RNA_seq" ,
      "assaymicrowell_seq" ,
      "assay10x_5_v2" ,
      "assay10x_5_v1" ,
      "assaySmart_seq2"
    )
  ) |>
  filter(parameter |> str_detect("^assay")) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  mutate(parameter = parameter |> str_remove("assay")) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>

  # Summarise and rank
  with_groups(parameter, ~ .x |> mutate(gran_mean = mean(v_effect))) |>
  arrange(desc(v_effect)) |>
  mutate(rank = formatC(1:n(), width = 3, format = "d", flag = "0")) |>

  mutate(cell_type_harmonised_label = glue("{rank}__{cell_type_harmonised}"))

source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7c78b50dce501fc7ce0b2a8d8efd3aded91134aa/color_cell_types.R")

# This figure shows the effects of the differential variability
plot_variability_error_bar_per_cell =
  res_relative_for_variability_plot |>

  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

  # Select the most variable cell types
  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |> with_groups(parameter, ~ .x |> arrange(desc(v_effect)) |>  slice_head(n=3)) |>

  ggplot(aes(v_effect, cell_type_harmonised_label)) +
  geom_errorbar(aes(xmin=v_lower, xmax=v_upper, color=cell_type_harmonised), ) +
  #geom_point() +
  facet_grid(fct_reorder(parameter, gran_mean) ~ ., scales = "free") +
  scale_color_manual(values = color_array) +
  scale_y_discrete(  labels = function(x) x |> str_remove("^[0-9]+__")   ) +
  xlab("Variability") +
  ylab("Cell group") +
  guides(color = "none") +
  theme_multipanel

# Get density plot for variability effects overall
plot_variability_density =
  res_relative_for_variability_plot |>

  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

	# Sample from effect distribution
  mutate(v_sd = (v_upper-v_effect)/qnorm(0.95)) |>
  mutate(distribution = map2(v_effect, v_sd, ~ rnorm(10, mean = .x, sd = .y))) |>
  unnest(distribution) |>
	
	# Plot
  ggplot(aes(distribution, color = parameter)) +
  geom_density(  alpha = 0.3) +
  stat_central_tendency(aes(color = parameter), type = "median", linetype = 2) +
  scale_color_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all("-", " ") |> str_remove("'")) |>  deframe() ) +
  guides(color = "none") +
  theme_multipanel


# Build plots with patchwork
plot =
  (
    ((
    	(
	    	(
	    		plot_assay_PCA |
	      wrap_heatmap(plot_heatmap_relative_assay, padding = unit(c(-30, 0, -3, -30), "points" ))
	    	) + plot_layout(  width = c(1,2) )
    	) /
    		plot_spacer()
    	) + plot_layout(  height = c(3.5,1) ) ) |
      ( ( plot_variability_density / plot_variability_error_bar_per_cell ) +  plot_layout(heights = c(1,4) ) )
  ) +
  plot_layout( guides = 'collect', widths = c(3, 1) ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")


ggsave(
	output_pdf,
  plot = plot,
  units = c("mm"),
  width = 183 ,
  height = 77 ,
  limitsize = FALSE
)

