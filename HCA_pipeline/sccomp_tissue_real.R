# Rscript /home/users/allstaff/mangiola.s//PostDoc/immuneHealthyBodyMap/HCA_pipeline/sccomp_tissue_real.R FALSE /home/users/allstaff/mangiola.s/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/input_common.rds /home/users/allstaff/mangiola.s/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_real_FALSE.rds /home/users/allstaff/mangiola.s/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_real_FALSE_blood.rds /home/users/allstaff/mangiola.s/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1/tissue_real_FALSE_proportion_adjusted.rds



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
input_blood_proportion = args[[4]]
output_file_2 = args[[5]]




res_relative =
	readRDS(input_file) |>
	
	mutate(is_lymphoid = tissue_harmonised %in% c(  	"blood",
																									 "lymph node",
																									 "bone",
																									 "spleen",
																									 "thymus")) |> 
	#--------------------#
	#Filter so samples are not experimentally enriched
	#----------------------#

	filter(file_id != "e756c34a-abe7-4822-9a35-55ef12270247") |>
	filter(file_id != "ca4a7d56-739b-4e3c-8ecd-28704914cc14") |>
	filter(
		file_id != "59dfc135-19c1-4380-a9e8-958908273756" |
			!tissue_harmonised %in% c("intestine small", "intestine large")
	) |>
	
	# Filter extremes
	nest(data = -c(.sample, tissue_harmonised, is_lymphoid)) |>
	mutate(n_immune = map_int(data, ~ .x |> filter(is_immune) |> nrow())) |> 
	mutate(n__NON_immune = map_int(data, ~ .x |> filter(!is_immune) |> nrow())) |> 
	
	# Filter samples which not include non immune cells
	filter(n_immune > 0) |> 
	
	# Filter samples which include > 90% of immune cells
	filter(is_lymphoid | (n_immune / (n_immune + n__NON_immune )) < 0.75) |> 
	
	unnest(data) |>
	
	#--------------------#
	# Filter low confidence
	#----------------------#
	filter(confidence_class %in% c(1, 2, 3)) |>
	
	# Fix groups
	unite("group", c(tissue_harmonised , file_id), remove = FALSE) |> 

  # Estimate
  sccomp_glm(
    formula_composition = ~ 0 + tissue_harmonised + sex + ethnicity_simplified  + age_days + assay_simplified + (tissue_harmonised | group),
    formula_variability = ~ 0 + tissue_harmonised + ethnicity_simplified + assay_simplified,
    .sample, cell_type_harmonised,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
  )

res_relative |> saveRDS(output_file_1)

# Remove unwanted variation
res_relative |>
		remove_unwanted_variation(~ 0 + tissue_harmonised, ~ 0 + tissue_harmonised) |>
		saveRDS(output_file_2)

