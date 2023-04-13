## Loading required packages:

library(dplyr)
library(dbplyr)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
library(scater)
library(scuttle)
library(Seurat)
library(Rtsne)
library(HCAquery) 
library(sccomp)

library(purrr) 
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(glue)

## Load metadata

meta.anno <- readRDS("../Data/metadata_annotated.rds")
meta.anno.sub <- meta.anno %>% filter(cell_type_harmonised!="non_immune") %>% filter(!is.na(confidence_class))
rm(meta.anno)

meta.anno.sub$confidence_class <- as.character(meta.anno.sub$confidence_class)
meta.anno.sub$.sample <- as.character(meta.anno.sub$.sample)


## glm for full model with age as random effect
estimate_age_days <-
  meta.anno.sub |>
  filter(!is.na(age_days)) |>
  unite("sample", c(.sample, assay, ethnicity, tissue_harmonised, sex, age_days), remove=FALSE) |>
  unite("re_age_days", c(dataset_id, age_days), remove = FALSE) |>
  sccomp_glm(
    formula_composition = ~ 0 + age_days + sex + tissue_harmonised + assay + ethnicity + 
      (age_days | re_age_days),
    formula_variability = ~ 1,
    sample, confidence_class,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    enable_loo = TRUE,
    prior_mean_variable_association = list(intercept = c(3.6539176, 2), slope = c(-0.5255242, 0.6), standard_deviation = c(20, 40))
  )

save(estimate_age_days, file="../Outputs/sccompestimate_age_days.Rdata")


# glm without age
estimate_age_days_0 <-
  meta.anno.sub |>
  filter(!is.na(age_days)) |>
  unite("sample", c(.sample, sex, tissue_harmonised, assay, ethnicity), remove=FALSE) |>
  sccomp_glm(
    formula_composition =  ~ 0 + sex + tissue_harmonised + assay + ethnicity,
    formula_variability = ~ 1,
    sample, confidence_class,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    enable_loo = TRUE,
    prior_mean_variable_association = list(intercept = c(3.6539176, 2), slope = c(-0.5255242, 0.6), standard_deviation = c(20, 40))
  )

save(estimate_age_days_0, file="../Outputs/sccompestimate_age_days_0.Rdata")




