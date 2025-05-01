library(Seurat)
library(ggplot2)

##
## - plot each features and counts for each count matrix 
## - using Seurat objects
## - helps to see which genes to filter due to abnormally high or low counts
##

directories <- list("filtered_feature_bc_matrix_a1",
                    "filtered_feature_bc_matrix_a2",
                    "filtered_feature_bc_matrix_b1",
                    "filtered_feature_bc_matrix_b2",
                    "filtered_feature_bc_matrix_c1",
                    "filtered_feature_bc_matrix_c2",
                    "filtered_feature_bc_matrix_d1",
                    "filtered_feature_bc_matrix_d2")

names <- list("condition_a_sample_1",
              "condition_a_sample_2",
              "condition_b_sample_1",
              "condition_b_sample_2",
              "condition_c_sample_1",
              "condition_c_sample_2",
              "condition_d_sample_1",
              "condition_d_sample_2")

i <- 1
for (dir in directories) {
  counts <- Read10X(dir)
  dim(counts)
  
  # create seurat object
  obj <- CreateSeuratObject(counts=counts)
  
  # get mito percent
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = '^MT-')
  
  # create plots
  # check for empty droplets and multiplets
  print(VlnPlot(obj, features = c("nFeature_RNA")) + 
                xlab("nFeature_RNA") + 
                ggtitle(names[i]) +
                scale_x_discrete(labels = NULL))
  print(VlnPlot(obj, features = c("nCount_RNA")) + 
                xlab("nCount_RNA") + 
                ggtitle(names[i]) + 
                scale_y_continuous(limits = c(0, 27000), breaks = seq(0, 27000, 5000)) +
                scale_x_discrete(labels = NULL))
  
  # check to make sure that counts and features are correlated
  print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(names[i]))
  
  # check to make sure high % mitochondrial transcripts are associated with low UMI counts
  print(FeatureScatter(obj, feature1 = "percent.mt", feature2 = "nCount_RNA") + ggtitle(names[i]))

  i = i + 1
}

# clear environment to conserve memory
rm(list=ls())