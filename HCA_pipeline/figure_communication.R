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

## from http://tr.im/hH5A

# get_metadata_local(cache_directory = "/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.1.sqlite") |>  filter(.cell == "GCACTAATCCTGGCTT_TSP1_endopancreas_1") |> select(.cell, file_id, lineage_1)


differential_composition_age_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_absolute_FALSE.rds"
data_for_immune_proportion_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/input_absolute.rds"
proportions_age_absolute_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/age_absolute_FALSE_proportion_adjusted.rds"

softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }
  
  exp(x - logsumexp(x))
}

dropLeadingZero <-
  function(l) {
    stringr::str_replace(l, '0(?=.)', '')
  }
S_sqrt <- function(x) {
  sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
  x ^ 2 * sign(x)
}
S_sqrt_trans <-
  function()
    scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

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

## from http://tr.im/hH5A
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)

data_for_immune_proportion_relative_file = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/input_relative.rds"
data_for_immune_proportion_relative = readRDS(data_for_immune_proportion_relative_file)

tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()


source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")

cell_type_color = 
  data_for_immune_proportion |> 
  pull(cell_type_harmonised) |> 
  unique() |> 
  get_cell_type_color()

names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")





# Communication
cdc_weight_fit = 
  readRDS("~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/cdc_gene_fit_random_intercept_file_resolved_cells_with_probability.rds") 

S_sqrt <- function(x) {
  sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
  x ^ 2 * sign(x)
}
S_sqrt_trans <-
  function()
    scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

lines_cell_type_communication = 
  cdc_weight_fit |>
  unnest(estimate_95) |> 
  filter((`l-95% CI` * `u-95% CI`)>0) |> 
  add_count(cell_type) |> 
  filter(n>3) |> 
  filter(!cell_type %in% c("non_immune", "immune_unclassified")) |> 
  with_groups(cell_type, ~ .x |> mutate(m = median(Estimate))) |> 
  mutate(cell_type = cell_type |> str_remove("phage") |> str_replace("th1/th17", "th1/17") |> str_replace("mono", "mn")) |> 
  
  # CLean names
  mutate(cell_type = cell_type |> 
           str_replace("tcm", "cm") |> 
           str_replace("cd4 th", "Th") |> 
           str_replace("memory", "mem") |> 
           str_replace("naive", "nv")
  ) |> 
  
  ggplot(aes(Estimate, fct_reorder(cell_type, m, .desc = TRUE) )) + 
  geom_vline(xintercept = 0, linetype="dashed", color="grey") + 
  geom_boxplot(outlier.shape = NA, width = 0) + theme(legend.position = "bottom", color="grey60") +
  geom_point(size = 0.2) +
  ylab("Cell type") +
  xlab("Fold changes for communication patterns") +
  scale_x_continuous(trans="S_sqrt") +
  theme_multipanel +
  theme(axis.text.y = element_text(angle=-40, hjust = 1, vjust = 1))


# Volcano plot
volcano_plot_communication =
  cdc_weight_fit  |>
  unnest(estimate_95) |> 
  filter(cell_type %in% c("cdc", "macrophage", "treg")) |> 
  
  mutate(significant = probability<0.05) |>
  mutate(DB = case_when(significant~DB)) |>
  mutate(gene = case_when(significant~gene)) |>
  ggplot(aes(Estimate, probability, label=gene)) +
  geom_point(aes(color=tissue_harmonised, size=significant)) +
  ggrepel::geom_text_repel(size = 1.5, max.overlaps = 20 ) +
  geom_hline(yintercept = 0.05, linetype="dashed", color="darkgrey", size = 0.2) +
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey", size = 0.2) +
  scale_size_discrete(range = c(0, 0.5)) +
  #scale_color_manual(values = c("#F1A340", "#998FC3"), na.value = "black") +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
  scale_x_continuous(trans="S_sqrt") +
  scale_color_manual(values = tissue_color) +
  facet_wrap(~cell_type, nrow=1) +
  xlab("Association communication strenght/age") +
  ylab("Probability H0") +
  theme_multipanel

draw_cellchat_circle_plot = function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
                                      targets.use = NULL, remove.isolate = FALSE, top = 1, top_absolute = NULL, weight.scale = T,
                                      vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15,
                                      vertex.label.cex = 0.8, vertex.label.color = "black", edge.weight.max = NULL,
                                      edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE,
                                      edge.label.color = "black", edge.label.cex = 0.8, edge.curved = 0.2,
                                      shape = "circle", layout = in_circle(), margin = 0.2, vertex.size = NULL,
                                      arrow.width = 1, arrow.size = 0.2)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  options(warn = -1)
  
  if(!is.null(top_absolute)) {
    thresh = top_absolute
    net[abs(net) < thresh] <- 0
  }
  
  thresh <- stats::quantile(as.numeric(net) %>% abs %>% .[.>0], probs = 1 - top)
  
  net[abs(net) < thresh] <- 0
  
  if(sum(net)==0) return(NULL)
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if(length(idx)>0){
      net <- net[-idx, ,drop=FALSE]
      net <- net[, -idx, drop=FALSE]
    }
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g),
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g),
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  
  if (is.null(color.use)) {
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$color <- color.use[igraph::V(g)]
  }
  else{
    igraph::V(g)$frame.color <- color.use[igraph::V(g) |> names()]
    igraph::V(g)$color <- color.use[igraph::V(g) |> names()]
  }
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(abs(igraph::E(g)$weight))
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + abs(igraph::E(g)$weight)/edge.weight.max *
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * abs(igraph::E(g)$weight)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  
  igraph::E(g)$color =
    circlize::colorRamp2(seq(max(abs(igraph::E(g)$weight)), -max(abs(igraph::E(g)$weight)), length.out =11), RColorBrewer::brewer.pal(11, "RdBu"))(igraph::E(g)$weight) %>%
    grDevices::adjustcolor(alpha.edge)
  
  
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[
      which(edge.start[, 2] == edge.start[,1])
    ] <- loop.angle[
      edge.start[
        which(edge.start[, 2] == edge.start[, 1]), 1
      ]
    ]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape,
       layout = coords_scale,
       margin = margin,
       #vertex.label.dist = label.dist,
       vertex.label.degree = label.locs,
       vertex.label.family = "Helvetica",
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 0.8)
  }
  gg <- recordPlot()
  return(gg)
}

pdf("~/PostDoc/CuratedAtlasQueryR/dev/circles_communication.pdf", width = 10, height = 33)
par(mfrow=c(10, 4))
plot_circles    =
  readRDS("~/PostDoc/CuratedAtlasQueryR/dev/cdc_gene_AND_cell_fit_random_intercept_file_resolved_cells.rds")  |>
  mutate(Estimate = Estimate/Est.Error/2) |>
  
  #filter(cell_from=="cdc") |>
  filter(cell_from!="non_immune" & cell_to != "non_immune") |>
  filter(cell_from!="immune_unclassified" & cell_to != "immune_unclassified") |>
  #filter((`l-95% CI` * `u-95% CI`)>0) |>
  
  nest(gene_data = -c(gene, DB, tissue_harmonised)) |>
  
  # Filter significant
  inner_join(
    cdc_weight_fit  |>
      unnest(estimate_95) |> 
      filter(cell_type %in% c("cdc", "macrophage", "treg")) |> 
      
      filter( probability<0.05) |> 
      distinct(gene, DB, tissue_harmonised)
  ) |> 
  left_join(
    CellChatDB.human$interaction|>
      #filter(pathway_name |> str_detect("CD23|MHC-I$|ICAM|ADGRE5|SELPLG|ITGB2")) |>
      nest(data = -pathway_name) |>
      mutate(ligand_receptor = map_chr(
        data,
        ~ glue("({.x$ligand |> unique() |> sort() |>  str_c(collapse = ',')}) -> \n ({.x$receptor |> unique() |> sort() |>  str_c(collapse = ',')})")
      )) |>
      mutate(ligand_receptor = ligand_receptor |>
               str_replace("HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G", "HLA-A-C/E-G") |>
               str_replace("CD8A,CD8B,CD8B2", "CD8A/B/B2") |>
               str_replace("ITGAL,ITGAL_ITGB2,ITGAM_ITGB2,ITGAX_ITGB2", "ITGAL/M/X_ITGB2") |>
               str_replace("RAET1E,RAET1F,RAET1G", "RAET1E-G") |>
               str_replace("CD94:NKG2A,CD94:NKG2C,CD94:NKG2E", "CD94:NKG2A/C/E") |>
               str_replace("KIR2DL1,KIR2DL2,KIR2DL3,KIR2DS1,KIR2DS4,KIR3DL1,KIR3DL2,KIR3DL3,KIR3DS1,KLRC1,KLRC2,KLRK1", "KIR2DL/DS/C") |>
               str_replace("LILRB1,LILRB2", "LILRB1/2") |>
               str_replace("ITGAM_ITGB2,ITGAV_ITGB3,ITGAX_ITGB2", "ITGAM/X_ITGB2,ITGAV_ITGB3")
      ),
    by=c("gene" = "pathway_name")
  ) |>
  mutate(title = glue("{tissue_harmonised}\n{ligand_receptor}")) |>
  # filter(gene %in% c("CD23", "MHC-I", "ICAM", "ADGRE5", "SELPLG", "ITGB2")) |>
  mutate(plot = map2(
    gene_data, title,
    ~ {
      
      cell_types = .x |> select(cell_from, cell_to) |> gather() |> pull(value) |>  unique()
      
      expand_grid(cell_from = cell_types, cell_to = cell_types) |>
        left_join(
          .x |>
            dplyr::select(cell_from, cell_to, Estimate)
        ) |>
        tidyr::complete(cell_from, cell_to, fill= list(Estimate = 0)) |>
        pivot_wider(names_from = cell_to, values_from =  Estimate) |>
        tidybulk:::as_matrix(rownames = cell_from) |>
        multiply_by(2) %>%
        
        # Avoid unknown error
        { 
          x = (.)
          x[x==0] = 0.000001
          x
        } |>
        
        draw_cellchat_circle_plot(
          edge.width.max = 4,
          remove.isolate = TRUE,
          #top = 0.2,
          arrow.width = 3,
          arrow.size = 0.3,
          edge.weight.max = 8,
          color.use = cell_type_color,
          title.name = .y
        )
    }
  ))
dev.off()



plot_communication = 
  lines_cell_type_communication + 
  volcano_plot_communication + 
  plot_spacer() + 
  plot_layout(guides = 'collect', widths = c(0.15, 0.55, 0.25)) &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )

ggsave(
  "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/figure_communication.pdf",
  plot = plot_communication,
  units = c("mm"),
  width = 183 ,
  height = 65 ,
  limitsize = FALSE
)
