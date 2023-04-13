start.time <- Sys.time()

# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(tidyverse)
library(stringr)
library(Seurat)
library(tidyseurat)

print("commandArgs")
args <- commandArgs(trailingOnly=TRUE)

print("Print args")
args


print("Get number_of_cells from args[1]")
number_of_cells <- as.numeric(gsub("--number_of_cells=", "", args[1]))


print("Get memory from args[2]")
memory <- gsub("--memory=", "", args[2])
memory <- as.numeric(gsub("G", "", memory))
memory <- memory/1024 ## Convert MB to GB

print("Get cores from args[3]")
cores <- gsub("--cores=", "", args[3])
cores <- as.numeric(cores)


print(paste("Number of cells =", number_of_cells))
print(paste("Memory requested (GB) =", memory))
print(paste("CPUs per task requested =", cores))

filename_part <- paste0("ncells", number_of_cells)



print("Read: variable_genes.rds")
variable_features <- readRDS("variable_features__CD4.rds")  
length(variable_features) 

# Path to the RDS files
file_directory = "Data/"

print("seurat_object")
seurat_object =
  dir(file_directory, full.names = TRUE)  |>
  enframe(value="file") |>
  mutate(seurat = map(file, ~ {            #purr::map(vector, function)
    seu =
      readRDS(.x)[variable_features,] |> 
      slice_sample(n=number_of_cells) |>
      select(-contains("name"))
    seu[["RNA"]] = NULL
    seu
  })) |>
  pull(seurat)


names(seurat_object) = dir(file_directory, full.names = TRUE)

print("Number of cells in each Batch")
lapply(seurat_object, function(x) ncol(x))

# seurat_object |> saveRDS("seurat_object.rds")
# #
# seurat_object = readRDS("seurat_object.rds")



# dataset using these features
print("PCA step")
seurat_object <- lapply(X = seurat_object, FUN = function(x) {
  x <- NormalizeData(x)
  x <- ScaleData(x, features = variable_features, verbose = FALSE)
  x <- RunPCA(x, features = variable_features, verbose = FALSE)
})

gc()

# To access the parallel version of functions in Seurat
print("future/sequential")

library(future)
plan(sequential)
plan("multicore", workers = cores)
options(future.globals.maxSize = memory * 1000 * 1024^2)


# WHERE THE PARALLEL - ancoring pairwise
print("immune.anchors")
immune.anchors <-
  
  seurat_object |>
  
  FindIntegrationAnchors(
    anchor.features = variable_features,
    reduction = "rpca"
  )

gc()
rm(seurat_object)

immune.anchors |> saveRDS(paste0("immune.anchors_", filename_part, ".rds"))


# immune.anchors = readRDS("immune.anchors_ncells800.rds")

# MAYBE  THE PARALLEL - intergation
print("immune.combined")
immune.combined <- IntegrateData(anchorset = immune.anchors, k.weight = 50)

immune.combined |> saveRDS(paste0("immune.combined_", filename_part, ".rds"))


end.time <- Sys.time()

elapsed.time <- difftime(time1=end.time, time2=start.time, units="secs")
print(elapsed.time)

sessionInfo()

