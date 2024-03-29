library(tidyverse)
library(forcats)
library(CuratedAtlasQueryR)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
library(ggforestplot)
library(ggrepel)
library(tidybulk)
library(ggforce)
library(ggpubr)
library(tidyHeatmap)
library(ComplexHeatmap)

# # Read arguments
# args = commandArgs(trailingOnly = TRUE)
# file_for_annotation_workflow = args[[2]]
# output_path = args[[3]]

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


# Load results
res_absolute = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_absolute_FALSE.rds")
data_for_immune_proportion = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/input_absolute.rds")
res_absolute_adjusted = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_absolute_FALSE_proportion_adjusted.rds")

# Define which are lymphoid organs
lymphoid_organs = c("blood", "spleen", "bone", "thymus", "lymph node")

# This function is used to format ggplots to save space
dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }

# Scale axis for ggplot
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)

#------------------------------#
# Analyses of immune cellularity proportion of immune cells in a tissue
#------------------------------#

# Remove unwanted variation from the immune cellularity
# Including effect of Age, sex, ethnicity, technology, 
# and random effects line datasets
# This is used to overlay boxplots on the observed proportions
res_generated_proportions =
	res_absolute |>
	sccomp_replicate(~ 0 + tissue_harmonised + sex + ethnicity  + age_days, number_of_draws = 100 ) |>
	filter(is_immune=="TRUE") |>
	left_join(
		data_for_immune_proportion |>
			dplyr::select(.sample, tissue_harmonised) |>
			distinct()
	) |>
	
	# Subset for each sample
	with_groups(tissue_harmonised, ~ .x |> sample_n(200, replace = T))
 

# Plot observed proportions of immune cells across non-lymphoid tissues
# And overlaying boxplots of the estimated proportion distribution
# Excluding the unwanted variation
plot_immune_proportion_dataset =
  data_for_immune_proportion |>
  
  # Stats
  dplyr::count(.sample, tissue_harmonised, is_immune, file_id) |>
  with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n), sum = sum(n))) |>
  filter(is_immune=="TRUE") |>
  with_groups(tissue_harmonised, ~ .x |> mutate( median_proportion = mean(proportion))) |>
  
  # Add multilevel proportion medians
  left_join(
  	res_generated_proportions |>
      with_groups(tissue_harmonised, ~ .x |> summarise(median_generated = median(generated_proportions, na.rm = TRUE)))
  ) |>
  
  # Label lymphoid organs
  mutate(is_lymphoid = tissue_harmonised %in% c("blood", "spleen", "bone", "thymus", "lymph node")) |>
  
  clean_names() |>
  
  # Arrange
  arrange(is_lymphoid, median_generated) %>%
  mutate(tissue_harmonised = factor(tissue_harmonised, levels = unique(.$tissue_harmonised))) |>
  
  # Cap
  mutate(sum = sum |> pmax(500)) |>
  mutate(sum = sum |> pmin(100000)) |>
  
  # Plot
  ggplot(aes( proportion, tissue_harmonised)) +
  ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_jitter(aes(size = sum, color=file_id), width = 0) +
  geom_boxplot(aes(generated_proportions, fct_reorder(tissue_harmonised, median_generated)), color="black", data =
  						 	res_generated_proportions |>
                 with_groups(tissue_harmonised, ~ .x |> mutate(median_generated = median(generated_proportions, na.rm = TRUE))) |>
                 clean_names(),
               fill = NA, outlier.shape = NA, lwd=0.2
  ) +
  guides(color="none") +
  scale_size(trans = "sqrt", range = c(0.1, 1.5), limits = c(500, 100000)) +
  scale_color_manual(values = dittoSeq::dittoColors()) +
  scale_x_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
  xlab("Immune proportion (sqrt)") +
  ylab("Tissue") +
  theme_multipanel


# Calculate the median composition for the abundance/variability plot
median_composition =
  res_absolute |>
  filter(factor == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) |>
  pull(c_effect) |>
  median()

# Calculate the median variability for the abundance/variability plot
median_variability =
  res_absolute |>
  filter(factor == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) |>
  pull(v_effect) |>
  median()

# Get estimates
coefficients_regression =
	res_absolute |>
	filter(factor == "tissue_harmonised") |>
	filter(c_effect<1) |>
	filter(is_immune == "TRUE") |>
	filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) %>%
	lm( v_effect ~ c_effect, data = .) %$%
	coefficients

# Scatter plot of abundance vs variability per tissue 
# Used for tissue landscape Figure panel B
res_for_plot =
  res_absolute |>
  filter(factor == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  mutate(tissue_harmonised = parameter |> str_remove("tissue_harmonised")) |>
  #mutate(intercept = coefficients_regression[1], slope = coefficients_regression[2]) |>
  

	# Define diversity in abundance and variability across tissues, far from the madian
	# Used for tissue landscape Figure panel B
  arrange(desc(c_effect)) |>
  mutate(	c_significant = ( row_number() <=7 | row_number() >= n() -3 ) & !tissue_harmonised %in%	c("blood", "lymph node", "spleen", "bone")	) |>
  arrange(desc(v_effect)) |>
  mutate(	v_significant = ( row_number() <=3 | row_number() >= n() -3) & !tissue_harmonised %in%	c("blood", "lymph node", "spleen", "bone")) |>
  
  # Define quadrants
  mutate(quadrant = case_when(
    c_effect > median_composition & v_effect > median_variability & (c_significant | v_significant) ~ "Hot and variable",
    c_effect > median_composition & v_effect < median_variability & (c_significant | v_significant) ~ "Hot and consistent",
    c_effect < median_composition & v_effect > median_variability & (c_significant | v_significant) ~ "Cold and variable",
    c_effect < median_composition & v_effect < median_variability & (c_significant | v_significant) ~ "Cold and consistent"
    
  )) |>
  
  # Limit the values
  mutate(
    v_effect = pmax(v_effect, -5),
    c_effect = pmin(c_effect, 1.5),
    v_lower = pmax(v_lower, -5),
    v_upper = pmax(v_upper, -5),
    c_lower = pmin(c_lower, 1.5),
    c_upper = pmin(c_upper, 1.5),
  )

# Plot abundance/variability plot
# Used for tissue landscape Figure panel B
plot_abundance_variability =
  res_for_plot |>
  ggplot(aes(c_effect, v_effect, label = parameter)) +
  geom_vline(xintercept = median_composition, linetype = "dashed",  alpha = 0.5, lwd = 0.2) +
  geom_hline(yintercept = median_variability, linetype = "dashed",  alpha = 0.5, lwd = 0.2) +
  geom_errorbar(aes(ymin = v_lower, ymax = v_upper, color = v_significant), alpha = 0.5, lwd = 0.2) +
  geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_significant), alpha = 0.5, lwd = 0.2) +
  geom_point(aes(fill = quadrant), shape = 21, stroke = 0.2) +
  #geom_text_repel(data =  res_for_plot |> filter(!v_significant & !c_significant)) +
  geom_text_repel(
    data = res_for_plot |>
      filter(v_significant | c_significant),
    size = 2
  ) +
  xlab("Composition") +
  ylab("Variability") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
  scale_fill_brewer(palette = "Set1") +
  theme_multipanel

# Clean environment
rm( res_absolute_adjusted, res_for_plot )
gc()



#------------------------------#
# Analyses of immune composition
#------------------------------#

# Read results
res_relative = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_relative_FALSE.rds")
data_for_immune_proportion_relative = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/input_relative.rds")
res_relative_adjusted = readRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_relative_FALSE_proportion_adjusted.rds")

# Get dataset for the tSNE 
# Figure tissue landscape panel C
observed_proportion_PCA_df =
	data_for_immune_proportion_relative |>
  
  # Mutate days
  filter(development_stage!="unknown") |>
  
  add_count(tissue_harmonised) |>
  filter(n > 5) |>
  dplyr::select(-n) |>
	
  dplyr::count(.sample, cell_type_harmonised, tissue_harmonised, assay, sex, file_id) |>
  with_groups(.sample, ~ .x |> mutate(observed_proportion = n/sum(n))) |>
  tidyr::complete(nesting(.sample, tissue_harmonised, assay, sex, file_id), cell_type_harmonised, fill = list(observed_proportion = 0)) |>
  reduce_dimensions(.sample, cell_type_harmonised, observed_proportion, method="tSNE", action="get")

# Plot for immune proportions WITH unwanted variation
# Coloured by tissue
observed_proportion_PCA_tissue =
  observed_proportion_PCA_df |>
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = tissue_harmonised), shape=21, stroke = NA, size=0.2) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none") +
  ggtitle("Tissue") +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Plot for immune proportions WITH unwanted variation
# Coloured by technology (technical variation)
observed_proportion_PCA_batch =
  observed_proportion_PCA_df |>
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = file_id), shape=21, stroke = NA, size=0.2) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none")  +
  ggtitle("Dataset") +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Plot for immune proportions with NO unwanted variation
# Including Age, ethnicity, sex, technology and random effects (datasets)
# Coloured by technology (technical variation)
adjusted_proportion_PCA =
	res_relative_adjusted |> 
  left_join(
    data_for_immune_proportion_relative |>
      distinct(.sample, tissue_harmonised, assay, file_id, sex, ethnicity)
  ) |>
  add_count(tissue_harmonised) |>
  filter(n > 5) |>
  reduce_dimensions(.sample , cell_type_harmonised, adjusted_proportion, method="UMAP", action="get") |>
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(fill = tissue_harmonised), shape=21, stroke = NA, size=0.2) +
  #ggdensity::geom_hdr_lines(aes(color = tissue_harmonised)) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  ggtitle("Adjusted tissue") +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Data for PCA plot where each point if a tissue
# Used for figure of tissue landscape panel D
 res_absolute_for_PCA =
  res_relative |>
  filter(factor == "tissue_harmonised") |>
  mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  rename(feature=cell_type_harmonised) |>
  bind_rows(
    res_absolute |>
      filter(factor=="tissue_harmonised") |>
      mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
      filter(is_immune=="TRUE") |>
      mutate(is_immune = "immune") |>
      rename(feature=is_immune)
  ) |>
  filter(parameter != "skeletal_muscle") |>
  reduce_dimensions(
    parameter, feature, c_effect,
    method="PCA", action="get", scale=FALSE,
    transform = identity
  )

# PCA plot for figure of tissue landscape panel D
plot_tissue_PCA =
  res_absolute_for_PCA |>
  ggplot(aes(PC1, PC2, label = parameter)) +
  geom_point(aes(fill = parameter), shape=21, stroke = 0.2) +
  ggrepel::geom_text_repel(size=2, max.overlaps = 10, min.segment.length = unit(10, "pt")) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none")  +
  theme_multipanel


tissue_color =
	data_for_immune_proportion_relative |>
	distinct(tissue_harmonised ) |>
	arrange(tissue_harmonised) |>
	mutate(color = dittoSeq::dittoColors()[1:n()]) |>
	deframe()

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")

# Color coding for tissue
cell_type_color = 
	data_for_immune_proportion |> 
	pull(cell_type_harmonised) |> 
	unique() |> 
	get_cell_type_color()
names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")

# Dataset used for the heatmap of immune composition across tissues
# These proportions are adjusted by unwanted variation
# Including Sex, ethnicity, age, technology, and random effects (datasets)
df_heatmap_relative_organ_cell_type =
  res_relative |>
  filter(factor == "tissue_harmonised") |>
  mutate(tissue_harmonised =   parameter ) |>
  clean_names() |>
  mutate(tissue =   tissue_harmonised ) |>
  mutate(cell_type = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  filter(cell_type !="immune_unclassified") |>
  
  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  ) |>
  
  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>
  
  with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |> mutate(top_3 = row_number() |> between(1, 3))) |>
  
  # Cell type abundance
  with_groups(tissue, ~ .x |> mutate(proportion = softmax(c_effect))) |>
  with_groups(cell_type, ~ .x |>  mutate(`Mean diff` = mean(proportion, na.rm = TRUE))) |>
  mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))

# heatmap of immune composition across tissues
# These proportions are adjusted by unwanted variation
# Including Sex, ethnicity, age, technology, and random effects (datasets)
plot_circle_relative_tissue =
  df_heatmap_relative_organ_cell_type |>
  filter(cell_type != "non_immune") |> 
	with_groups(parameter, ~ .x |> mutate(proportion = softmax(c_effect))) |> 
  mutate(proportion_label = proportion |> round(3) |> dropLeadingZero())  |> 
  	
	# Use tidyHeatmap
  heatmap(
    tissue, cell_type, proportion,
    palette_value = circlize::colorRamp2(seq(-8, -2, length=5), viridis::magma(5)),
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    show_heatmap_legend = FALSE,
    row_km = 2, 
    transform = car::logit
  ) |>
  split_rows(6) |>
  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  
  annotation_tile(
    tissue, show_legend = FALSE,
    palette =
      df_heatmap_relative_organ_cell_type |>
      distinct(tissue, tissue_color) |>
      arrange(tissue) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    cell_type, show_legend = FALSE,
    palette =
      df_heatmap_relative_organ_cell_type |>
      distinct(cell_type, cell_type_color)  |>
      arrange(cell_type) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    count_tissue, show_legend = FALSE,
    size = unit(0.2, "cm"),
    palette = circlize::colorRamp2(c(0, 5, 10, 15), viridis::magma(4))
  )  |>
	layer_text(.value = proportion_label, .size = 4)


# Compose plot with patchwork
first_line_first_column =
  
  (plot_immune_proportion_dataset|plot_tissue_PCA) +
	plot_layout( guides = 'collect', width = c(10,6) ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")



p =
  (
  	first_line_first_column /
    wrap_heatmap(plot_circle_relative_tissue, padding = unit(c(-67, -10, -0, -30), "points" ))
  ) + 
	plot_layout( guides = 'collect', heights = c(6,9) ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")




ggsave(
  "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/figure_tissue.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 230 ,
  limitsize = FALSE
)


plot_circle_relative_tissue |> 
	save_pdf("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/figure_tissue_heatmap.pdf", width = 183, height=110, units="mm")


