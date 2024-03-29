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

differential_composition_ethnicity_absolute =
  readRDS(input_file) |>

  # Drop only-immune organs
  filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone")) |>
  mutate(is_immune = as.character(is_immune)) |>

  filter(development_stage!="unknown") |>

  # Scale days
  mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>

	# Create one-to-many grouping for multilevel modelling
  unite("tissue_harmonised_ethnicity", c(tissue_harmonised , ethnicity_simplified), remove = FALSE) |>

  # Estimate
  sccomp_glm(
    formula_composition = ~ 0 + ethnicity_simplified + tissue_harmonised + sex  + age_days +  assay_simplified  + (ethnicity_simplified | tissue_harmonised_ethnicity),
    formula_variability = ~ 0 + ethnicity_simplified + tissue_harmonised + sex,
    .sample, is_immune,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
  )

differential_composition_ethnicity_absolute |>
  saveRDS(output_file_1)



# Counts RUV Absolute
differential_composition_ethnicity_absolute |>
  remove_unwanted_variation(~ 0 + ethnicity_simplified + (ethnicity_simplified | tissue_harmonised_ethnicity), ~ ethnicity_simplified) |>
  saveRDS(output_file_2)

