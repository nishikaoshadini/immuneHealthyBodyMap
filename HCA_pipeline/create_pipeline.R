library(tidyverse)
library(glue)
library(here)
library(stringr)
library(Seurat)
library(tidyseurat)

args = commandArgs(trailingOnly=TRUE)
run_directory = args[[1]]

vast_directory = "/vast/projects/RCP/human_cell_atlas"
R_code_directory = glue("/home/users/allstaff/mangiola.s//PostDoc/immuneHealthyBodyMap/HCA_pipeline")
root = "/home/users/allstaff/mangiola.s//PostDoc/immuneHealthyBodyMap"
tab = "\t"
metadata_sqlite = glue("{vast_directory}/metadata_annotated_0.1.5.sqlite")

harmonised_annotation = glue("{root}/cell_metadata_with_harmonised_annotation.rds")
lineage_df = "~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv"

run_directory |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)


commands = c()

# Create input
commands =
  commands |>
  c(
    glue("CATEGORY=create_input\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
    glue("{run_directory}/input_common.rds {run_directory}/input_absolute.rds {run_directory}/input_relative.rds:{harmonised_annotation} {metadata_sqlite}\n{tab}Rscript {R_code_directory}/create_common_data.R {metadata_sqlite} {lineage_df} {run_directory}/input_common.rds {run_directory}/input_absolute.rds {run_directory}/input_relative.rds")
  )

# Estimate blood contamination
commands =
  commands |>
  c(
    glue("CATEGORY=blood_contamination\nMEMORY=30024\nCORES=10"),
    glue("{run_directory}/blood_fit.rds {run_directory}/blood_contamination.rds:{run_directory}/input_relative.rds {run_directory}/input_absolute.rds\n{tab}Rscript {R_code_directory}/estimate_blood_contamination.R {run_directory}/input_relative.rds {run_directory}/input_absolute.rds {run_directory}/blood_fit.rds {run_directory}/blood_contamination.rds")
  )

# estimate
estimate_commands =
	tibble(
		factor = c("tissue", "sex", "ethnicity", "assay", "age")
	) |>
	expand_grid(
		modality = c("absolute", "relative"),
		filter_blood = c("FALSE") #c("TRUE", "FALSE")
	) |>
	filter(!(factor == "assay" & modality=="absolute")) |>
	mutate(analysis = glue("{factor}_{modality}_{filter_blood}")) |>
	mutate(r_script = glue("sccomp_{factor}_{modality}.R")) |>
	mutate(input = glue("input_{modality}.rds")) |>
	rowwise() |>
	mutate(command = if_else(
		filter_blood == "TRUE",
		glue("{run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds:{run_directory}/{input}\n{tab}Rscript {R_code_directory}/{r_script} {filter_blood} {run_directory}/{input} {run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds"),
		glue("{run_directory}/{analysis}.rds {run_directory}/{analysis}_proportion_adjusted.rds:{run_directory}/{input}\n{tab}Rscript {R_code_directory}/{r_script} {filter_blood} {run_directory}/{input} {run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds")
	)
	)

commands =
  commands |>
  c(
    glue("CATEGORY=estimate\nMEMORY=60024\nCORES=20"),
    estimate_commands |>
    	pull(command)
  )

# Plots
commands =
	commands |>
	c(
		glue("CATEGORY=plots\nMEMORY=30024\nCORES=10"),
		glue("{run_directory}/assay.rds {run_directory}/assay.pdf:{run_directory}/assay_relative_FALSE.rds {harmonised_annotation} {metadata_sqlite}\n{tab}Rscript {R_code_directory}/figure_assay.R {run_directory}/assay_relative_FALSE.rds {harmonised_annotation} {metadata_sqlite} {run_directory}/assay.rds {run_directory}/assay.pdf")

	)

commands |>
  write_lines(glue("{run_directory}/pipeline.makeflow"))
