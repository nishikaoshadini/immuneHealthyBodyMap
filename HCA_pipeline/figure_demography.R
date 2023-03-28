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
library(scales)
library(ggplot2)


result_directory = "~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1"


# Calculate softmax from an array of reals
softmax <- function (x) {
	logsumexp <- function (x) {
		y = max(x)
		y + log(sum(exp(x - y)))
	}
	
	exp(x - logsumexp(x))
}

# This function is used to format ggplots to save space
dropLeadingZero <-
	function(l) {
		stringr::str_replace(l, '0(?=.)', '')
	}

# Scale axis for ggplot
S_sqrt <- function(x) {
	sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
	x ^ 2 * sign(x)
}
S_sqrt_trans <-
	function()
		scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

# This function shorten names of cell types for visualisation purposes
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


# Read inout files 
data_for_immune_proportion_absolute_file = glue("{result_directory}/input_absolute.rds")
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)
data_for_immune_proportion_relative_file = glue("{result_directory}/input_relative.rds")
data_for_immune_proportion_relative = readRDS(data_for_immune_proportion_relative_file)

# Color coding for tissue
tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")
cell_type_color = 
  data_for_immune_proportion |> 
  pull(cell_type_harmonised) |> 
  unique() |> 
  get_cell_type_color()
names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")


#------------------------------#
# Sex analyses for immune cellularity
#------------------------------#

# Read results
differential_composition_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE.rds")
differential_composition_sex_relative_file = glue("{result_directory}/sex_relative_FALSE.rds")
proportions_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE_proportion_adjusted.rds")
differential_composition_sex_relative = readRDS(differential_composition_sex_relative_file)  
differential_composition_sex_absolute = readRDS(differential_composition_sex_absolute_file)
proportions_sex_absolute_adjusted = readRDS(proportions_sex_absolute_file)

# Get parameter draws for relative abundance 
# to manually plot the uncertainty
draws_abundance = 
  readRDS(differential_composition_sex_absolute_file)   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE")

# Get parameter draws for variability 
# to manually plot the uncertainty
draws_variability = 
  readRDS(differential_composition_sex_absolute_file)   |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE") |> filter(parameter=="sexmale")

# Plot effect of composition and variability
# For immune cellularity (proportion of immune cells)
# Wtih uncertainty
plot_sex_absolute_1D =
  tibble(
    Variability = draws_variability |> pull(.value),
    Abundance = draws_abundance |> pull(.value)
  ) |>
  tidybulk::as_matrix() |>
  bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
  coord_flip() +
  xlab("Effect male immune cellularity") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=90, hjust = 0.5))

# Create dataset to create the mannequin heatmap of the 
# Tissues with differential immune cellularity
sex_absolute_organ_tissue =
  differential_composition_sex_absolute |> 
  test_contrasts(
    contrasts =
      differential_composition_sex_absolute |>
      filter(parameter |> str_detect("___sex")) |>
      distinct(parameter) |>
      mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      filter(contrast |> str_detect("_female", negate = TRUE)) |> 
      deframe( ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  filter(is_immune == "TRUE")

# Draw the color palette for the mannequin heatmap of the 
# Tissues with differential immune cellularity
colors_palette_for_organ_abundance =
  sex_absolute_organ_tissue |>
  select(parameter, c_effect) |> 
  mutate(color = circlize::colorRamp2(
    seq(1.45,-1.45, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )(c_effect)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>
  pull(color) |>
  scales::show_col(	cex_label = 0.5	)

# Draw the boxplot for significant changes of the 
# Tissues with differential immune cellularity
plot_sex_absolute_organ_boxoplot_adjusted =
  proportions_sex_absolute_adjusted |>
  left_join(
    data_for_immune_proportion |>
      distinct(.sample, tissue_harmonised, ethnicity, sex,tissue, file_id)
  ) |>
  inner_join(
    sex_absolute_organ_tissue |> 
      filter(c_FDR<0.07) |> 
      separate(parameter, c("tissue_harmonised", "sex"), sep="_") |>  
      distinct(tissue_harmonised)
  ) |> 
  
  mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
  filter(is_immune =="TRUE") |>
  ggplot(aes(sex, adjusted_proportion )) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
  geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
  facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
  guides(color = "none") +
  ylab("Adjusted proportion") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))



#------------------------------#
# Sex analyses for immune composition
#------------------------------#

# Volcano plot of cell type difference overall between sexes
volcano_relative_sex = 
  differential_composition_sex_relative |>
  
  test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
  filter(cell_type_harmonised != "non_immune") |> 
  
  filter(parameter == "sexmale") |> 
  mutate(naive_experienced = case_when(
    cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
    cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
  )) |> 
  mutate(significant = c_FDR<0.05) |> 
  mutate(cell_type_harmonised = case_when(c_FDR<0.05~cell_type_harmonised)) |> 
  ggplot(aes(c_effect, c_FDR)) + 
  
  geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
  
  geom_point(aes(color=naive_experienced, size=significant)) +
  ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
  scale_x_continuous(trans="S_sqrt") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_size_discrete(range = c(0, 0.5)) +
  guides(size="none", color="none") +
  theme_multipanel

#------------------------------#
# Ethnicity analyses for immune cellularity
#------------------------------#


# Read results

differential_composition_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE.rds")
proportions_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE_proportion_adjusted.rds")
differential_composition_ethnicity_relative_file = glue("{result_directory}/ethnicity_relative_FALSE.rds")
differential_composition_ethnicity_relative = readRDS(differential_composition_ethnicity_relative_file)
differential_composition_ethnicity_absolute = readRDS(differential_composition_ethnicity_absolute_file)
proportions_ethnicity_tissue_absolute_adjusted = readRDS(proportions_ethnicity_absolute_file)

 # Test differential immune cellularity 
# using contrasts for each ethnicity compared with the 
# Average of the other three
data_for_ethinicity_absolute_plot =
  differential_composition_ethnicity_absolute |>
  test_contrasts(
    c(
      African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
      Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
      European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
      Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
    )
  ) |>
  
  filter(parameter %in% c("African",  "Hispanic", "European", "Chinese")) |>
  filter(is_immune == "TRUE")

# Plot effect uncertainty for abundance
draws_abundance = 
  differential_composition_ethnicity_absolute   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(
    African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
    Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
    European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
    Asian = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
  )) |> 
  filter(is_immune=="TRUE")

# Plot effect uncertainty for varability
draws_variability = 
  differential_composition_ethnicity_absolute  |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(
    African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
    Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
    European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
    Asian = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
  )) |> 
  filter(is_immune=="TRUE") |> 
  filter(parameter %in% c("Hispanic", "African", "European", "Asian")) 

# Plot the effects for abundance and variability 
plot_ethnicity_absolute_1D =
  left_join(
    draws_variability |> select(.value, parameter, .draw) |> rename(Variability = .value),
    draws_abundance |> select(.value, parameter, .draw) |> rename(Abundance = .value),
  ) |>
  select(-.draw) |> 
  select(parameter, everything()) |> 
  nest(data = -parameter) |> 
  mutate(plot = map2(
    data, parameter,
    ~ .x |> 
      tidybulk::as_matrix() |>
      bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
      coord_flip() +
      xlab("Effect male immune cellularity") +
      xlim(-0.75, 0.75) +
      theme_multipanel +
      theme(axis.text.x = element_text(angle=90, hjust = 0.5)) +
      ggtitle(.y)
  )) |> 
  rownames_to_column() |> 
  mutate(plot = map2(
    plot, rowname, 
    ~ if(.y>1) .x + 
      theme(
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()
      ) else .x)) |> 
  pull(plot) |> 
  wrap_plots(nrow=1)
  
# Test the difference in immune cellularity across ethnicities
# Per tissue
ethnicity_absolute_organ =
  differential_composition_ethnicity_absolute |>
  
  # Find stats of random effect with groups
  test_contrasts(
  	
  	# Create contrasts dinamically
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
      select(tissue_harmonised, contrast) |> 
      separate(tissue_harmonised, c("tissue_harmonised", "ethnicity"), sep="_") |> 
      add_count(tissue_harmonised) |> mutate(contrast = glue("({contrast})")) |>  
      pivot_wider(names_from = ethnicity, values_from = contrast) |> 
      
      # Filter where there are at least 3 ethnicities
      mutate(n_minus_1 = n-1) |> 
      filter(n_minus_1>1) |>
      
      # Build contrasts
      rowwise() |> 
      mutate(
        other_Chinese = c(European, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        other_European = c(Chinese, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        other_African = c(European, Chinese, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
        `other_Hispanic or Latin American` = c(European, African, Chinese) |> str_subset(".+") |> str_c(collapse=" + "),
      ) |> 
      mutate(
        contrast_Chinese = glue("1/{n_minus_1}*({other_Chinese}) - {Chinese}"),
        contrast_European = glue("1/{n_minus_1}*({other_European}) - {European}"),
        contrast_African = glue("1/{n_minus_1}*({other_African}) - {African}"),
        contrast_Hispanic= glue("1/{n_minus_1}*({`other_Hispanic or Latin American`}) - {`Hispanic or Latin American`}"),
      ) |> 
      rename(
        this_Chinese = Chinese,
        this_European = European,
        this_African = African,
        this_Hispanic = `Hispanic or Latin American`
      ) |> 
      select(tissue_harmonised, starts_with(c("contrast_", "this_"))) |>
      pivot_longer(cols = starts_with(c("contrast_", "this_")), names_to = c("type", "ethnicity"), names_sep = "_") |> 
      pivot_wider(names_from = type, values_from = value) |>
      filter(!is.na(this)) |> 
      select(tissue_harmonised, ethnicity, contrast) |> 
      unite("tissue_harmonised", c(tissue_harmonised, ethnicity)) |> 
      deframe( )
  )



# Palette for the maniquin heatmap
# About the tissue-level difference across ethnicities
# For immune cellularity
ethnicity_absolute_organ |>
  filter(is_immune == "TRUE") |>
  select(parameter, c_effect, c_FDR) |> 
  separate(parameter, c("tissue", "ethnicity"), remove = FALSE, sep="_") |> 
  mutate(color = circlize::colorRamp2(
    seq(2,-2, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )(c_effect)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>
  
  # group
  nest(data = -ethnicity) |> 
  mutate(dummy = map(data, ~ scales::show_col(.x$color,	cex_label = 0.5	))) |> 
  pull(data)

# Plot the boxplot visualisation for the plot above
plot_ethnicity_absolute_organ_boxoplot_adjusted =
  proportions_ethnicity_tissue_absolute_adjusted |>
  left_join(
    data_for_immune_proportion |>
      distinct(.sample, tissue_harmonised, ethnicity, tissue, file_id)
  ) |>
  inner_join(
    ethnicity_absolute_organ |> 
      filter(c_FDR<0.05) |>  
      separate(parameter, c("tissue", "ethnicity"), remove = FALSE, sep="_") |> 
      distinct(tissue) |> rename(tissue_harmonised = tissue)
  ) |>
  mutate(ethnicity = ethnicity |> str_remove(" or Latin American") ) |>
  mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
  filter(is_immune =="TRUE") |>
  ggplot(aes(ethnicity, adjusted_proportion )) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
  geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
  facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
  guides(color = "none") +
  ylab("Adjusted proportion") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))


#------------------------------#
# Ethnicity analyses for immune composition
#------------------------------#

# Test difference in immune composition at the body
# Comparing each ethnicity with the average of the others
data_for_ethinicity_relative_plot =
  differential_composition_ethnicity_relative |>
  test_contrasts(
    c(
      African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
      Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
      European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
      Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
    ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  
  filter(parameter %in% c("African",  "Hispanic", "European", "Chinese"))

# Plot volcano of the differences at the immune composition level
# (for each cell type)
volcano_relative_ethnicity = 
  data_for_ethinicity_relative_plot |>
  filter(cell_type_harmonised != "non_immune") |> 
  mutate(naive_experienced = case_when(
    cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
    cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
  )) |>
  mutate(significant = c_FDR<0.05) |> 
  mutate(cell_type_harmonised = case_when(c_FDR<1e-2~cell_type_harmonised)) |> 
  ggplot(aes(c_effect, c_FDR)) + 
  geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
  geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
  
  geom_point(aes(color=naive_experienced, size=significant)) +
  ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
  facet_wrap(~parameter, nrow=1) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
  scale_x_continuous(trans="S_sqrt") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_size_discrete(range = c(0, 0.5)) +
  theme_multipanel


# Clean environment
rm(differential_composition_ethnicity_relative , differential_composition_ethnicity_absolute )
gc()

# Compose plot with patchwork
plot = 
  ((
    ((
      plot_sex_absolute_1D |
        plot_spacer() |
        plot_sex_absolute_organ_boxoplot_adjusted |
        volcano_relative_sex
    ) + plot_layout(width = c(12, 22, 12, 30))) /
      
      plot_spacer() /
      
      (plot_ethnicity_absolute_1D  + plot_layout(width = c(2,2,2,2)) ) /
      
      volcano_relative_ethnicity 
      
  ) + plot_layout( height = c(55, 18, 20, 50))
) |
  ((
    plot_spacer() / plot_ethnicity_absolute_organ_boxoplot_adjusted / plot_spacer()
  )  + plot_layout( height = c(85, 50, 38)) ) +
  plot_layout(width = c(84, 97)) &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )




ggsave(
  glue("~/PostDoc/CuratedAtlasQueryR/dev/sccomp_on_HCA_0.2.1/figure_demography.pdf"),
  plot = plot,
  units = c("mm"),
  width = 183 ,
  height = 180 ,
  limitsize = FALSE
)
