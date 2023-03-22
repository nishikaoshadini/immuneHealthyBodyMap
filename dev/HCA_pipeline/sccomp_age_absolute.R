library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")


args = commandArgs(trailingOnly=TRUE)
filter_blood = args[[1]]
input_file = args[[2]]
output_file_1 = args[[3]]
output_file_1_blood = args[[4]]
output_file_2 = args[[5]]

immune_non_immune_differential_composition_age =
  readRDS(input_file) |>

		# Drop only-immune organs
		# filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone", "thymus")) |>
		mutate(is_immune = as.character(is_immune)) |>

		# Mutate days
  filter(development_stage!="unknown") |>
		unite("group", c(tissue_harmonised , file_id), remove = FALSE) |>

	# Filter embrios
	filter(age_days_original > 365) |>

		sccomp_glm(
			formula_composition = ~age_days + tissue_harmonised + sex + ethnicity + assay + (1 | group) + (age_days | tissue_harmonised),
			formula_variability = ~age_days + tissue_harmonised ,
			.sample, is_immune,
			check_outliers = F,
			approximate_posterior_inference = FALSE,
			cores = 20,
			mcmc_seed = 42,
			verbose = T,
			prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
		)

immune_non_immune_differential_composition_age |>

		saveRDS(output_file_1)

immune_non_immune_differential_composition_age |>
    remove_unwanted_variation( ~ age_days, ~ age_days) |>
    saveRDS(output_file_2)
