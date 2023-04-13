# Loading required packages:
  

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
library(corrplot)
library(reshape2)
library(viridis)
library(tidyverse)
library(loo)


## Load the metadata

meta.anno <- readRDS("../Data/metadata_annotated.rds")

## Filter immune cells

meta.anno.sub <- meta.anno %>% filter(cell_type_harmonised!="non_immune") %>% filter(!is.na(confidence_class))
rm(meta.anno)


## Distribution of confidence class across across ethnicity (Plot E)

library(dittoSeq)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

df.ethncty <- meta.anno.sub |> count(ethnicity, confidence_class)
colnames(df.ethncty) <- c("ethnicity", "confidence_class", "count")
df.ethncty <- df.ethncty |> group_by(ethnicity) |> mutate(prop=count/sum(count))

df.ethncty <- df.ethncty |>
  nest(data = -ethnicity) |>
  mutate(rank_quantity = map_dbl(data, ~ .x |> filter(confidence_class ==1) |> pull(prop))) |>
  unnest(data) 

#Re-factor ethnicity
df.ethncty$ethnicity2 <- factor(df.ethncty$ethnicity, levels = c("African", "African American", "African American or Afro-Caribbean",
                                                                 "Asian", "Chinese", "European", "Han Chinese", "Hispanic or Latin American", 
                                                                 "unknown"))
levels(df.ethncty$ethnicity2) <- c("African", "AfroAmerican", "AfroAmerican/Caribbean",
                                   "Asian", "Chinese", "European", "HanChinese", "Hispanic/LatAmerican", 
                                   "unknown")


plot.ethnicity <- ggplot(df.ethncty) + aes(x=fct_reorder(ethnicity2,rank_quantity, .desc = TRUE), y=prop, fill=confidence_class) + geom_col() +
  scale_fill_viridis_c(direction=-1) + theme_multipanel +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  labs(fill='confidence_class') + ylab("Proportion") + xlab("Ethnicity") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15))

ggsave("../Plots/plot.ethnicity.pdf",
       units = c("mm"),
       width = 55 ,
       height = 56, #max=230mm
       limitsize = FALSE)




## Distribution of confidence class across across assay (Plot F)

df.assay <- meta.anno.sub |> count(assay, confidence_class)
colnames(df.assay) <- c("assay", "confidence_class", "count")
df.assay <- df.assay |> group_by(assay) |> mutate(prop=count/sum(count))

df.assay <- df.assay |>
  nest(data = -assay) |>
  mutate(rank_quantity = map_dbl(data, ~ .x |> filter(confidence_class ==1) |> pull(prop))) |>
  unnest(data) 

# Re-factor assay
df.assay$assay2 <- factor(df.assay$assay, levels = c("10x 3' transcription profiling", "10x 3' v1", "10x 3' v2",
                                                     "10x 3' v3", "10x 5' v1", "10x 5' v2", "10x scATAC-seq",
                                                     "microwell-seq", "sci-RNA-seq", "scRNA-seq", "Seq-Well", "Slide-seq",
                                                     "Smart-seq2", "Visium Spatial Gene Expression"))
levels(df.assay$assay2) <- c("10x 3'", "10x 3' v1", "10x 3' v2",
                             "10x 3' v3", "10x 5' v1", "10x 5' v2", "10xscATAC-seq",
                             "microwell-seq", "sci-RNA-seq", "scRNA-seq", "Seq-Well", "Slide-seq",
                             "Smart-seq2", "Visium")


plot.assay <- ggplot(df.assay) + aes(x=fct_reorder(assay2, rank_quantity, .desc = TRUE), y=prop, fill=confidence_class) + geom_col() +
  scale_fill_viridis_c(direction=-1) + theme_multipanel +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  labs(fill='confidence_class') + ylab("Proportion") + xlab("Assay")

ggsave("../Plots/plot.assay.pdf",
       units = c("mm"),
       width = 75 ,
       height = 57,
       limitsize = FALSE)



## Distribution of confidence class across across assay cell type (Plot C)

df.celltype2 <- meta.anno.sub |> count(cell_type_harmonised, confidence_class)
colnames(df.celltype2) <- c("cell_type_harmonised", "confidence_class", "count")
df.celltype2 <- df.celltype2 |> group_by(cell_type_harmonised) |> mutate(prop=count/sum(count))

# Since cell type "mast" doesn't have confidence_class = 1, 
# we add a new row to it with count = 0 for confidence_class = 1
# in order to order the cell_types based on values of confidence_class 1
# in the following code.

tmp <- data.frame("mast", 1, 0, 0)
names(tmp) <- colnames(df.celltype2)

df.celltype2 <- rbind(df.celltype2, tmp)


df.celltype2 <- df.celltype2 |>
  nest(data = -cell_type_harmonised) |>
  mutate(rank_quantity = map_dbl(data, ~ .x |> filter(confidence_class ==1) |> pull(prop))) |>
  unnest(data) 

plot.celltype2 <- ggplot(df.celltype2) + aes(x=fct_reorder(cell_type_harmonised, rank_quantity, .desc = TRUE), y=prop, fill=confidence_class) + geom_col() +
  scale_fill_viridis_c(direction=-1) + theme_multipanel +
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  labs(fill='confidence_class') + ylab("Proportion") + xlab("Cell type")

ggsave("../Plots/plot.celltype2.pdf",
       units = c("mm"),
       width = 88 ,
       height = 48,
       limitsize = FALSE)




## Co-occurence of cell types across algorithms (Plot G)

tmp <- meta.anno.sub |> filter(confidence_class != 1)

#Obtain every possible cell type pair (136 pairs)
cellpair.mat <- combn(unique(tmp$cell_type_harmonised),2)

#Frequency of each cell type pair appears in rows of the 3 variables, 
#[cell_type_harmonised, cell_annotation_azimuth_l2 and cell_annotation_blueprint_singler]. 
tmp2 <- as.matrix(tmp[,c("cell_type_harmonised" ,"cell_annotation_azimuth_l2", "cell_annotation_blueprint_singler")])
out <- rep(0, ncol(cellpair.mat))
for (i in 1:ncol(cellpair.mat)) {
  tmp3 <- rep(0, nrow(tmp2))
  for (j in 1:nrow(tmp2)) {
    tmp3[j] <- all(cellpair.mat[,i] %in% tmp2[j,])
  }
  out[i] <- sum(tmp3)
}
out <- data.frame(t(cellpair.mat), out)
colnames(out) <- c("cell_type_1", "cell_type_2", "freq")
#save(out, file="../Outputs/freq_celltype_pair.Rdata")
load(file="../Outputs/freq_celltype_pair.Rdata")

#Order the frequency of cell type pairs in as descending order
out2 <- out[order(-out$freq), ]
out2$pair <- paste0(out2$cell_type_1,", ",out2$cell_type_2)

#Keep the frequencies > 300
df.cellpair <- out2[out2$freq>300,]

#Assign each cell type (columns: cell_type_1, cell_type_2) into major cell type groups:
#cd4 tcm, cd4 tem, cd4 naive, cd8 tem, cd8 tcm, cd8 naive, treg : Tcell
#cd14 mono, cd16 mono, macrophage, cdc : monocitic 
#b memory, b naive : Bcell
#pdc : pdc
#stem : stem
#mast : mast
#nk : nk
df.cellpair$cell_type_1_group <- df.cellpair$cell_type_1
df.cellpair$cell_type_1_group[df.cellpair$cell_type_1 %in% 
                                c("cd4 tcm", "cd4 tem", "cd4 naive", "cd8 tem", "cd8 tcm", "cd8 naive", "treg")] <- "Tcell"
df.cellpair$cell_type_1_group[df.cellpair$cell_type_1 %in% 
                                c("cd14 mono", "cd16 mono", "macrophage", "cdc")] <- "monocitic"
df.cellpair$cell_type_1_group[df.cellpair$cell_type_1 %in% c("b memory", "b naive")] <- "Bcell"

df.cellpair$cell_type_2_group <- df.cellpair$cell_type_2
df.cellpair$cell_type_2_group[df.cellpair$cell_type_2 %in% 
                                c("cd4 tcm", "cd4 tem", "cd4 naive", "cd8 tem", "cd8 tcm", "cd8 naive", "treg")] <- "Tcell"
df.cellpair$cell_type_2_group[df.cellpair$cell_type_2 %in% 
                                c("cd14 mono", "cd16 mono", "macrophage", "cdc")] <- "monocitic"
df.cellpair$cell_type_2_group[df.cellpair$cell_type_2 %in% c("b memory", "b naive")] <- "Bcell"

#Identify cell pairs with common and different major cell types:
df.cellpair$colour <- c(immunecell.col)
df.cellpair$colour[df.cellpair$cell_type_1_group != df.cellpair$cell_type_2_group] <- "indianred4"

df.cellpair$lineage <- "common lineage"
df.cellpair$lineage[df.cellpair$colour == "indianred4"] <- "different lineage"


#Plot
br <- c(0,2.5,5,7.5,10)
options(scipen=5)

plot.confoundcell <- ggplot(df.cellpair, aes(y=freq^0.2, x = reorder(pair, -freq), fill=lineage)) + 
  scale_fill_manual(values = c("#528B8B", "indianred4")) +
  geom_bar(stat="identity") + theme_multipanel + 
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  xlab("Cell type pair") + ylab("Counts") +
  scale_y_continuous(breaks = br, labels = round(exp(log(br)/0.2))) +
  guides(fill=guide_legend(ncol=1), shape = guide_legend(override.aes = list(size=0.2)))


ggsave("../Plots/plot.confoundcell.pdf",
       units = c("mm"),
       width = 140 ,
       height = 75,
       limitsize = FALSE)




## Autocorrelation between cell types

tmp <- meta.anno.sub |> filter(confidence_class != 4)
out <- table(tmp$.sample, tmp$cell_type_harmonised)
out <- data.frame(rbind(out))

#Data frame of the distribution (proportions) of samples per celltype
out.pr <- out/rowSums(out)

#Scale proportions to get rid of any zeros
compress_zero_one = function(y){ 
  n = length(y)
  (y * (n-1) + 0.5) / n
}

for (i in 1:ncol(out.pr)) {
  out.pr[,i] <- compress_zero_one(out.pr[,i])
}

#Calculate pearson correlation on logged data
logout <- log(out.pr)
cormat <- cor(logout, method="pearson")

get_lower_tri<-function(cormat){
  cormat[base::upper.tri(cormat)] <- NA
  return(cormat)
}

melted_cormat <- melt(get_lower_tri(cormat), na.rm = TRUE)

plot.corr <- ggplot(melted_cormat, aes(Var1, Var2)) +
  xlab('Cell type') +
  ylab('Cell type') +
  geom_tile(aes(fill = value), color='white') +
  scale_fill_gradient2(limits = c(-1,1), low = "steelblue", mid = "white", high = "darkred") +
  theme_multipanel +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee')) 

ggsave("../Plots/plot.corr.pdf",
       units = c("mm"),
       width = 60 ,
       height = 70,
       limitsize = FALSE)



#Calculate pearson correlation on non-logged data
cormat2 <- cor(out.pr, method="pearson")

get_upper_tri<-function(cormat){
  cormat[base::lower.tri(cormat)] <- NA
  return(cormat)
}

melted_cormat2 <- melt(get_upper_tri(cormat2), na.rm = TRUE)

plot.corr2 <- ggplot(melted_cormat2, aes(Var1, Var2)) +
  xlab('Cell type') +
  ylab('Cell type') +
  geom_tile(aes(fill = value), color='white') +
  scale_fill_gradient2(limits = c(-1,1), low = "steelblue", mid = "white", high = "darkred") +
  theme_multipanel +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee')) 

ggsave("../Plots/plot.corr2.pdf",
       units = c("mm"),
       width = 60 ,
       height = 70,
       limitsize = FALSE)



## Statistical modelling (Plot H)


# Ethnicity

load(file = "../Outputs/sccompestimate_ethnicity.Rdata")
load(file = "../Outputs/sccompestimate_ethnicity_0.Rdata")

# Compare models
comp.ethnicity <- loo_compare(
  estimate_ethnicity |> attr("fit") |> loo(),
  estimate_ethnicity_0 |> attr("fit") |> loo()
)

r.ethnicity <- comp.ethnicity[2,1]/comp.ethnicity[2,2]


# Tissue

load(file = "../Outputs/sccompestimate_tissue.Rdata")
load(file = "../Outputs/sccompestimate_tissue_0.Rdata")

comp.tissue <- loo_compare(
  estimate_tissue |> attr("fit") |> loo(),
  estimate_tissue_0 |> attr("fit") |> loo()
)

r.tissue <- comp.tissue[2,1]/comp.tissue[2,2]

# Assay

# When fitting a model for assay, use the file metadata_annotated_0.2.rds 
# which doesn't have assay in sample names
load(file = "../Outputs/sccompestimate_assay.Rdata")
load(file = "../Outputs/sccompestimate_assay_0.Rdata")

comp.assay <- loo_compare(
  estimate_assay |> attr("fit") |> loo(),
  estimate_assay_0 |> attr("fit") |> loo()
)

r.assay <- comp.assay[2,1]/comp.assay[2,2]


# Sex

load(file = "../Outputs/sccompestimate_sex.Rdata")
load(file = "../Outputs/sccompestimate_sex_0.Rdata")

comp.sex <- loo_compare(
  estimate_sex |> attr("fit") |> loo(),
  estimate_sex_0 |> attr("fit") |> loo()
)

r.sex <- comp.sex[2,1]/comp.sex[2,2]


# Age

load(file = "../Outputs/sccompestimate_age_days.Rdata")
load(file = "../Outputs/sccompestimate_age_days_0.Rdata")

comp.age <- loo_compare(
  estimate_age_days |> attr("fit") |> loo(),
  estimate_age_days_0 |> attr("fit") |> loo()
)

r.age <- comp.age[2,1]/comp.age[2,2]


# Plot

df.ratio <- data.frame(Variable = c("Sex", "Age", "Ethnicity", "Tissue", "Assay"), 
                       Ratio = -c(r.sex, r.age, r.ethnicity, r.tissue, r.assay))
#save(df.ratio, file="../Outputs/df.ratio.Rdata") 
load(file="../Outputs/df.ratio.Rdata") 

plot.ratio <- ggplot(df.ratio, aes(y=Ratio, x = Variable)) + 
  geom_bar(stat="identity", fill="grey") + theme_multipanel + 
  theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) +
  xlab("Variable") + ylab("-(elpd_diff / se_diff)") +
  geom_hline(yintercept=5, size=0.2, linetype=2) + ylim(c(0,30)) 

ggsave("../Plots/plot.ratio.pdf",
       units = c("mm"),
       width = 33 ,
       height = 42,
       limitsize = FALSE)
