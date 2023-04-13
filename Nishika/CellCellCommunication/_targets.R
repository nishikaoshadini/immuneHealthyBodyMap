# Created by use_targets().

# Load packages required to define the pipeline:
library(targets)

# Load custom functions and small input objects into the R session:
source("R/functions.R")

# Declare the packages that the targets need to run,
# as well as other settings such as the default storage format:
tar_option_set(packages = c("dplyr", "CuratedAtlasQueryR", "SingleCellSignalR", "Seurat", "stringr"),
               resources = tar_resources(clustermq = tar_resources_clustermq(template = list(num_cores = 48, memory = 512000)))
)# memory(MB), walltime(min)

# tar_make_clustermq() configuration
options(clustermq.scheduler = "slurm", clustermq.template = "clustermq.tmpl")

# Write the pipeline:
list(
  tar_target(name=metadata, command=split_metadata(), iteration = "group"),
  tar_target(name=signal, command=ccc_sample(metadata),  pattern=map(metadata))
)
