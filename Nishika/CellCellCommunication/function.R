# split_metadata() function keeps only the observations corresponding to confidence_class 1, 2 and 3. 
# It then filtered the cell types where the number of cells in each cell type is > 10 and 
# the cell types where the number of cells in each cell type is > 10 Within each sample.
# Finally, it subsets the filterd metadata into small groups based on the sample id (sample_).

split_metadata <- function(){
  
  my_df = 
    get_metadata() |>
    select( sample_, cell_type_harmonised) |>
    as_tibble()
  
  # List of the cell types where the number of cells in each cell type is > 10
  cell_type_with_many_cells = 
    my_df |>
    count(cell_type_harmonised) |>
    filter(n>10) |>
    pull(cell_type_harmonised)
  
  # List of the cell types where the number of cells in each cell type is > 10 Within each sample.
  cell_type_plus_sample_with_many_cells = 
    my_df |>
    count(sample_, cell_type_harmonised) |>
    filter(n>10)
  
  rm(my_df)
  gc()
  
  # Keep confidence_class 1, 2 and 3
  get_metadata() |>
    filter(confidence_class < 4) |> 
    
    # Filter
   filter(cell_type_harmonised %in% cell_type_with_many_cells) |>
    inner_join(cell_type_plus_sample_with_many_cells, copy = TRUE) |>
    
    collect() |> ## Retrieves data into a local tibble
    group_by(sample_) |>
    tar_group()
}


# ccc_sample() function inputs metadata per sample id. 
# It first obtain the seurat object for that sample id, normalizes it 
# and then calculates cell cell interaction scores and save the results.

ccc_sample <- function(metadata){
  
  ## Get seurat object per sample 
  srtsample <- metadata |> get_seurat(cache_directory = "/vast/projects/cellxgene_curated")
  
  ## Normalize data
  norm_srtsample <- NormalizeData(srtsample)
  
  ## Generate cellular interaction list per sample
  all.genes <- rownames(srtsample)
  data <- data.frame(srtsample[["originalexp"]]@data)
  cluster <- srtsample@meta.data[,"cell_type_harmonised"]
  cluster <- factor(cluster)
  cluster <- as.numeric(unclass(cluster))
  
  tar_dir(signal <- cell_signaling(data=data, genes=all.genes, cluster=cluster))
  #assign(paste0("signal", sampleid), cell_signaling(data=data, genes=all.genes, cluster=cluster))
  

  ## Save the cellular interaction list
  sampleid <- srtsample$sample_ |> unique()
  saveRDS(signal, file = paste0("RDSfiles/", sampleid, ".rds"))
  #save(list=paste0("signal", sampleid), file = paste0("RDSfiles/", sampleid, ".rds"))
}

