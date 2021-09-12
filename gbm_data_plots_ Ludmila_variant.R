# for unified and processed RNA-seq data
library(recount3)
# to normalize the RNA-seq data 
library(DESeq2) 
# for access to TCGA data
library(TCGAbiolinks)
# to look at the data
library(tidyverse)
# to visualize the mutation data
library(maftools)
# to create heatmaps
library(ComplexHeatmap)
scale2 <- function(mat, ...) {
  t(scale(t(mat), ...))
}

rse_gene <- 
  create_rse(
    subset(
      available_projects(),
      project == "GBM" & project_type == "data_sources"
    )
  )

assayNames(rse_gene)

assay(rse_gene, "counts") <- 
  transform_counts(rse_gene)

sample_sheet <-
  colData(rse_gene) %>%
  data.frame() %>%
  rownames_to_column("sample_id")


normalized_counts <- 
  vst(assays(rse_gene)$counts)

#View(assays(rse_gene)$counts)

#heatmaps
row_var <-
  rowVars(normalized_counts)

ht <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100))

ht


sample_sheet %>% 
  select(tcga.tcga_barcode, tcga.cgc_sample_sample_type) %>% 
  mutate(patient_id = str_extract(tcga.tcga_barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>% 
  group_by(patient_id) %>% 
  summarise(count = n(), 
            sample_type = paste(unique(sort(tcga.cgc_sample_sample_type)),
                                collapse = ", ")) %>% 
  filter(count > 1)



ha <-
  sample_sheet %>% 
  select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)


ht2 <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100), top_annotation = ha)
ht2




#Somatic mutations

library(TCGAbiolinks)
#data downloaded from  https://drive.google.com/file/d/1a_4mSFyDGQP_hXGGc6HHzluELv6C1gZr/view
maf <- readRDS('/home/marichka/Downloads/data')

plotmafSummary(maf = maf, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE, log_scale = TRUE)


tcga_mutation_burden <- 
  tcgaCompare(maf = maf, cohortName = "GBM")

oncoplot(maf = maf, top = 10)

lollipopPlot(maf, "EGFR", labelPos = 'all')
lollipopPlot(maf, "PTEN")
somaticInteractions(maf, top = 15, pvalue = c(0.01, 0.05))

#clinical data
sample_sheet %>%
  select(starts_with("tcga.gdc_cases.demographic."))

tcga_subtype_data <-
  TCGAquery_subtype(tumor = "gbm")

tcga_subtype_data %>%
  select(ends_with("subtype"))


library(ggplot2)
library(RColorBrewer)
library(ggprism)
library(patchwork)
library(magrittr)
library(rstatix)


df_p_val <- rstatix::t_test(tcga_subtype_data,  Age..years.at.diagnosis. ~ Gender) %>% 
  rstatix::add_x_position()

ggplot(data = tcga_subtype_data[complete.cases(tcga_subtype_data$Gender),], 
       aes(x = Gender, y = Age..years.at.diagnosis.)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "Set2") +
  xlab('Gender') + ylab('Age of diagnostics') + add_pvalue(df_p_val, y.position = 75)

colnames(tcga_subtype_data)
summary(tcga_subtype_data)
