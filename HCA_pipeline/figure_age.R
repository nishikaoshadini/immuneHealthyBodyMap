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
library(scales)
library(ggplot2)

## from http://tr.im/hH5A

# get_metadata_local(cache_directory = "/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.1.sqlite") |>  filter(.cell == "GCACTAATCCTGGCTT_TSP1_endopancreas_1") |> select(.cell, file_id, lineage_1)

# Set up files names 
differential_composition_age_absolute_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/age_absolute_FALSE.rds"
data_for_immune_proportion_absolute_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/input_absolute.rds"
proportions_age_absolute_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/age_absolute_FALSE_proportion_adjusted.rds"

# Calculate softmax from an array of reals
softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }

  exp(x - logsumexp(x))
}

# This function is used to format ggplots to save space
dropLeadingZero <-
	function(l) {
		stringr::str_replace(l, '0(?=.)', '')
	}

# Scale axis for ggplot
S_sqrt <- function(x) {
	sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
	x ^ 2 * sign(x)
}
S_sqrt_trans <-
	function()
		scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

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

# Read files
## from http://tr.im/hH5A
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)
data_for_immune_proportion_relative_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/input_relative.rds"
data_for_immune_proportion_relative = readRDS(data_for_immune_proportion_relative_file)

# Color coding for tissue
tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")
cell_type_color = 
  data_for_immune_proportion |> 
  pull(cell_type_harmonised) |> 
  unique() |> 
  get_cell_type_color()
names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")



# Threshold that equates to a linear increase of 1% from 20% to 21%
# This convoluted threshold gives a lay meaning to a increase in the softmax space 
# which is not linear and hard to grasp the meanin of a lay audience
# We use this threshold to test for significant effect bigger that it
FDR_threshold_1_percent_change_at_20_percent_baseline = 0.017
  
#------------------------------#
# Analyses of immune cellularity proportion of immune cells in a tissue
#------------------------------#

# Track of immune system in life
differential_composition_age = readRDS(differential_composition_age_absolute_file)

# Function that calculate the approximate proportional change from an effect in the iverse logit  space
# Again this is helpful for result interpretation
print_estimate_plus_minus = function(fit, contrasts_baseline, contrasts, contrasts_uncertainty){
	
	# 1.73 is the scaled 84
	
	baseline = 
		fit |> 
		test_contrasts(
			contrasts = contrasts_baseline,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		) |> 
		filter(is_immune == "TRUE") |> 
		pull(c_effect) 
	
	baseline = softmax(c(baseline, -baseline))[1]
	
	fit_contrasts = 
		fit |> 
		test_contrasts(
			contrasts = contrasts,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		) |> 
		filter(is_immune == "TRUE") 
	
	# Increase with credible interval
	estimate = 
		fit_contrasts |> 
		select(parameter, c_lower, c_effect, c_upper) |> 
		mutate(proportion_effect = map_dbl(c_effect, ~ softmax(c(.x, -.x))[1] - baseline )) |> 
		mutate(proportion_effect_10_years = proportion_effect/7.3) |> 
		select(parameter , proportion_effect_10_years) |> 
		arrange(parameter) |> 
		deframe()
	
	# Incertainty
	proportional_uncertainty = 
		fit |> 
		test_contrasts(
			contrasts = contrasts_uncertainty,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		)  |> 
		filter(is_immune == "TRUE") |> 
		mutate(uncertainty_in_percentage = c_upper/c_effect) |> 
		select(parameter , uncertainty_in_percentage) |> 
		arrange(parameter) |> 
		deframe()
	
	print(glue("{names(estimate)} {estimate} +/- {(estimate*proportional_uncertainty) - estimate}"))
}

# Age absolute stats in terms of proportion
# 1.73 is the equivalent of 84 years that jave been scaled to standard deviation of 1
differential_composition_age |> 
	print_estimate_plus_minus(
		contrasts_baseline = c(year_0 = "`(Intercept)`"),
		contrasts = c(year_84 = "`(Intercept)` + 1.73 * age_days"), 
		contrasts_uncertainty = c("age_days")
	)

# Plot change in immune cellularity through ageing
# generate points for the line trend in the plot
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

# Plot change in immune cellularity through ageing
# generate points for the line trend in the plot
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

# Plot change in immune cellularity through ageing
# generate points for the line trend in the plot
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


# Read proportion adjusted for unwanted variation, including sex, ethnicity, technology, and random effects
proportions_age_absolute = readRDS(proportions_age_absolute_file)

# Define life stages for result interpretation
life_stages = tibble(
  start = c(0, 2, 5, 13, 20, 40, 60),
  end = c(1, 4, 12, 19, 39, 59,100)
) |>
  mutate(stage = c("Infant", "Toddler", "Child", "Teen", "Adult", "Middle age", "Senior"))


# Add rectangles to the plot for life stages
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

# Figure 3B
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



# Plot association of immune cellularity with age, per tissue 
age_absolute_organ =
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
    test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
  ) |>
  filter(is_immune == "TRUE")

# Print statistics for percentage increase to be used in the paper
differential_composition_age |> 

	print_estimate_plus_minus(
		contrasts_baseline = c(baseline = "`(Intercept)`"),
		
		contrasts =
			differential_composition_age |>
			filter(parameter |> str_detect("___age_days")) |>
			distinct(parameter) |>
			mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
			mutate(contrast = glue("`(Intercept)` + (1.73 * ({contrast}))")) |> 
			mutate(tissue_harmonised = glue("year_84__{tissue_harmonised}")) |> 
			deframe(), 
		
			contrasts_uncertainty = 	differential_composition_age |>
															filter(parameter |> str_detect("___age_days")) |>
															distinct(parameter) |>
															mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
															tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+")  |> 
															deframe()
		)

# Plot the fold gain for SUPPLEMENTARY figures
adjusted_counts_for_cellularity_tissue_effects = 
	differential_composition_age |>
	remove_unwanted_variation( ~ age_days + tissue_harmonised + (age_days | tissue_harmonised), ~ age_days) |> 
	inner_join(data_for_immune_proportion |>
						 	tidybulk::pivot_sample(.sample) 
	)

adjusted_counts_for_cellularity_tissue_effects|> 
	filter(tissue_harmonised |> str_detect("brain|prostate|kidney")) |> 
	filter(is_immune == "TRUE") |> 
	ggplot(aes(age_days_original/365, adjusted_proportion)) + geom_point() + 
	facet_wrap(~tissue_harmonised) + scale_y_continuous(trans = "logit") +
	geom_smooth(method = "lm") + geom_hline(yintercept = 0.05) + geom_hline(yintercept = 0.15)

# Plot the fold gain decrease SUPPLEMENTARY figures
adjusted_counts_for_cellularity_tissue_effects  |> 
	filter(tissue_harmonised |> str_detect("heart|skin|adipose")) |> 
	filter(is_immune == "TRUE") |> 
	ggplot(aes(age_days_original/365, adjusted_proportion)) + geom_point() + 
	facet_wrap(~tissue_harmonised) + scale_y_continuous(trans = "logit") +
	geom_smooth(method = "lm", formula = y ~ x + I(x^2)) + geom_hline(yintercept = 0.05) + geom_hline(yintercept = 0.15)

# Color gradient based on effect to be used in the homan silhuette in biorender
# Figure 3
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


# Print significance global statistics to be used in the paper
count_significance_age_immune_load =
	age_absolute_organ_cell_type |>
	filter(is_immune=="TRUE") |>
	filter(c_FDR<0.05)


#------------------------------#
# Analyses of immune composition
#------------------------------#

# Load data
differential_composition_age_relative_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/age_relative_FALSE.rds"
proportions_age_relative_file = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/age_relative_FALSE_proportion_adjusted.rds"
differential_composition_age_relative =  readRDS(differential_composition_age_relative_file)
proportions_age_relative = readRDS(proportions_age_relative_file)

# Function that calculate the approximate proportional change from an effect in the inverse logit  space
# Again this is helpful for result interpretation
print_estimate_plus_minus_relative = function(fit, contrasts_baseline, contrasts, contrasts_uncertainty){
	
	# 1.73 is the scaled 84
	
	baseline = 
		fit |> 
		test_contrasts(
			contrasts = contrasts_baseline,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		) |> 
		select(cell_type_harmonised, parameter,  baseline = c_effect) |> 
		mutate(baseline = softmax(c(baseline)))


	fit_contrasts = 
		fit |> 
		test_contrasts(
			contrasts = contrasts,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		) 
	
	# Increase with credible interval
	estimate = 
		fit_contrasts |> 
		select(cell_type_harmonised, parameter, c_lower, c_effect, c_upper) |> 
		left_join(baseline) |> 
		with_groups(parameter, ~ .x |> mutate(proportion_effect = softmax(c_effect) - baseline)) |> 

		mutate(proportion_effect_10_years = proportion_effect/7.3) |> 
		
		mutate(fold_proportion_change = proportion_effect/baseline) |> 
		
		select(cell_type_harmonised, parameter , fold_proportion_change) |> 
		tidyr::unite("parameter", c(cell_type_harmonised, parameter)) |> 
		arrange(parameter) |> 
		deframe()
	
	# Incertainty
	proportional_uncertainty = 
		fit |> 
		test_contrasts(
			contrasts = contrasts_uncertainty,
			test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
		)  |> 
		mutate(uncertainty_in_percentage = c_upper/c_effect) |> 
		select(cell_type_harmonised, parameter , uncertainty_in_percentage) |> 
		tidyr::unite("parameter", c(cell_type_harmonised, parameter)) |> 
		arrange(parameter) |> 
		deframe()
	
	print(glue("{names(estimate)} {estimate} +/- {(estimate*proportional_uncertainty) - estimate}"))
}

# Print stats of global changes for immune composition in association with age 
differential_composition_age_relative |> 
	print_estimate_plus_minus_relative(
		contrasts_baseline = c(year_0 = "`(Intercept)`"),
		contrasts = c(year_84 = "`(Intercept)` + (1.73 * age_days)"), 
		contrasts_uncertainty = c("age_days")
	)


# Scale the x axis 
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

# Volcano of the global compositional changes
volcano_relative = 
 differential_composition_age_relative |>
	test_contrasts(test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline) |> 
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


# Print the fold change statistics to be used in the paper
differential_composition_age |> 
	print_estimate_plus_minus(
		contrasts_baseline = c(year_0 = "`(Intercept)`"),
		contrasts = c(year_84 = "`(Intercept)` + 1.73 * age_days"), 
		contrasts_uncertainty = c("age_days")
	)


# Get trend line to be used in the scatter plot of global compositional changes, plot_age_relative
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


# Scatter plot of the significant cell types which compositon change therough ageing
# The proportions are adjusted to exclude other effects including
# Sex, ethnicity, random effects (e.g. datasets) and technnology
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


# Fold change for cell type by organ
# These statistics are used in the paper result section
differential_composition_age_relative |> 
	
	print_estimate_plus_minus_relative(
		contrasts_baseline = 
			differential_composition_age_relative |>
			filter(parameter |> str_detect("___age_days")) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
			mutate(parameter = glue("__{parameter}")) |> 
			distinct(parameter, tissue_harmonised) |>
			mutate(contrast = glue("`(Intercept)` + `tissue_harmonised{tissue_harmonised}`") |> as.character()) |>
			select(parameter, contrast) |> 
			filter(parameter |> str_detect("adipose", negate = TRUE)) |> 
			deframe(),
		
		contrasts =
			differential_composition_age_relative |>
			filter(parameter |> str_detect("___age_days")) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
			distinct(parameter, tissue_harmonised) |>
			mutate(contrast = glue("age_days  + `{parameter}`") |> as.character()) |>
			mutate(contrast = glue("`(Intercept)` + `tissue_harmonised{tissue_harmonised}` + (1.73 * ({contrast}))")) |> 
			mutate(parameter = glue("__{parameter}")) |> 
			select(parameter, contrast) |> 
			filter(parameter |> str_detect("adipose", negate = TRUE)) |> 
			
			
			deframe(), 
		
		contrasts_uncertainty = 	
			differential_composition_age_relative |>
			filter(parameter |> str_detect("___age_days")) |>
			distinct(parameter) |>
			mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE)  |> 
			filter(parameter |> str_detect("adipose", negate = TRUE)) |> 
			mutate(parameter = glue("__{parameter}")) |> 
			
			deframe()
	)

# Prepare the dataset for drawing the heatmap of changes by tissues and celltypes, Figure 3
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
    test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
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

# Plot heatmap with tidyHeatmap, Figure 3
plot_heatmap_age_relative_organ_cell_type =

  df_heatmap_age_relative_organ_cell_type |>
  
  # Heatmap
  heatmap(
    tissue, cell_type, Difference,
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

# Save heatmap separately
plot_heatmap_age_relative_organ_cell_type |>
  save_pdf(
    filename = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/plot_heatmap_age_relative_organ_cell_type.pdf",
    width = 80*1.5, height = 60*1.5, units = "mm"
  )

# Color palette for the absolute effect sum for tissues
# This is used for colouring the human maniquine drawn with biorender
# And represents the most affected tissue through ageing for immune compositional changes
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


# Print the statistics for the maniquine above, to be used for the image of biorender
count_significance_age_cell_type =
  differential_composition_age_relative |>
  test_contrasts(test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline) |>
  filter(parameter=="age_days") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)

count_significance_age_cell_type_tissue =
  df_heatmap_age_relative_organ_cell_type |>
  count(c_FDR<0.05)

# Clean environment
rm(differential_composition_age_relative , differential_composition_age )
gc()

# Plot the association of cell residency with age
plot_trends_residency =
  readRDS("~/PostDoc/immuneHealthyBodyMap/residency_data.rds") |>
  
  filter(cell_type_harmonised |> str_detect("naive", negate = T)) |>
  separate(cell_type_harmonised, "cell_type_harmonised", sep = " ") |>
  
  # Filter significant
  inner_join(
    readRDS("~/PostDoc/immuneHealthyBodyMap/residency_estimates_random_effects.rds") |> 
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



# Compose plots with patchwork
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
  "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/figure_age.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 200 ,
  limitsize = FALSE
)



