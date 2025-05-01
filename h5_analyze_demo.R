library(Seurat)
library(RegenOrNoRegen)
library(garnett)
library(stringr)
library(tidyverse)
library(scater)
library(ggplot2)


setwd("/Users/emmagomez/code/axonogenesis/r_sim_new_data/h5_data/data")

# market matrix directories of sparse UMI count matrices
DIRECTORIES <- list("filtered_feature_bc_matrix_SARMKO_SHAM_5231_gex_molecule_info",
                    "filtered_feature_bc_matrix_SARMKO_SHAM_5234_gex_molecule_info",
                    "filtered_feature_bc_matrix_SARMKO_STROKE_5232_gex_molecule_info",
                    "filtered_feature_bc_matrix_SARMKO_STROKE_5233_gex_molecule_info",
                    "filtered_feature_bc_matrix_WT_SHAM_1_gex_molecule_info",
                    "filtered_feature_bc_matrix_WT_SHAM_2_gex_molecule_info",
                    "filtered_feature_bc_matrix_WT_STROKE_4_gex_molecule_info",
                    "filtered_feature_bc_matrix_WT_STROKE_5_gex_molecule_info")

NAMES <- list("SARMKO_SHAM_5231",
              "SARMKO_SHAM_5234",
              "SARMKO_STROKE_5232",
              "SARMKO_STROKE_5233",
              "WT_SHAM_1",
              "WT_SHAM_2",
              "WT_STROKE_4",
              "WT_STROKE_5")

# filtering global variables
N_FEATURE_LOWER_BOUND <- 500
N_COUNT_RNA_LOWER_BOUND <- 1000
N_COUNT_RNA_UPPER_BOUND <- 20000
N_FEATURE_UPPER_BOUNDS <- list()

# assign upper bounds, WT_SHAM_2 needs a slightly higher upper bound
for (dir in DIRECTORIES) {
  N_FEATURE_UPPER_BOUNDS[[dir]] <- 6000
}
N_FEATURE_UPPER_BOUNDS[["filtered_feature_bc_matrix_WT_SHAM_2_gex_molecule_info"]] <- 7000

i <- 1
# run the classification pipeline for each count matrix
for (dir in DIRECTORIES) {
  counts <- Read10X(dir)
  
  # create seurat object
  obj <- CreateSeuratObject(counts=counts)
  
  obj <- subset(obj,
                subset = nFeature_RNA > N_FEATURE_LOWER_BOUND & 
                         nFeature_RNA < N_FEATURE_UPPER_BOUNDS[[dir]] &
                         nCount_RNA > N_COUNT_RNA_LOWER_BOUND & 
                         nCount_RNA < N_COUNT_RNA_UPPER_BOUND
                )
  
  # normalize, find variable features, and scale
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  
  # run PCA dim reduction
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims=1:18)
  
  # classify cells
  print(SeuratLoad(obj, GeneIDType = "ENSEMBL"))
  obj <- SeuratReturn(obj, GeneIDType = "ENSEMBL")
  
  # get regeneration index
  regen_index <- table(obj$Regeneration_Index)
  regenerating <- unname(regen_index['Regenerating'])
  print(regen_index)
  
  # proportion df
  proportion_df <- as.data.frame(prop.table(regen_index))
  proportion_df$Class <- rownames(proportion_df)
  rownames(proportion_df) <- NULL
  colnames(proportion_df) <- c("Class", "Proportion")
  print(proportion_df)
  
  # bar chart
  chart_title <- paste("Proportions of Regenerating Cells:", NAMES[[i]])
  
  print(ggplot(proportion_df, aes(x = Class, y = Proportion, fill = Class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::percent(Proportion, accuracy = 0.01)), vjust =-0.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.001)) +
    ggtitle(chart_title) +
    theme_minimal() +
    ylab("Proportion") +
    xlab("Regeneration Class") +
    scale_fill_manual(values = c("Unknown" = "lightcoral", "NonRegenerating" = "#00bb38", "Regenerating" = "#629dff")))
  
  i = i + 1
}