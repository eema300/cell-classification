library(R.utils)
library(rhdf5)
library(data.table)
library(Matrix)
library(Seurat)

gzfile <- function(file) {
  R.utils::gzip(file, overwrite = TRUE)
}

files <- list("condition_a_sample_1.h5", "condition_a_sample_2.h5", 
              "condition_b_sample_1.h5", "condition_b_sample_2.h5",
              "condition_c_sample_1.h5", "condition_c_sample_2.h5",
              "condition_d_sample_1.h5", "condition_d_sample_2.h5")


for (file in files) {
  # peek at the file structure
  h5ls(file)
  
  # save relevant data
  # (barcode indices, feature indices, umis)
  barcode_idx <- h5read(file, "/barcode_idx")
  feature_idx <- h5read(file, "/feature_idx")
  umi <- h5read(file, "/umi")
  
  # (actual cell barcods, gene names, gene ids)
  barcodes <- h5read(file, "/barcodes")
  genes <- h5read(file, "/features/name")
  gene_ids <- h5read(file, "/features/id")
  
  feature_type <- rep("Gene Expression", length(gene_ids))
  
  # create map of gene-cell-umi indices
  dt <- data.table(
    gene = feature_idx + 1,
    cell = barcode_idx + 1,
    umi = umi
  )
  
  # check to see if any mitochondrial transcripts
  print(sum(grepl("^MT-", gene_ids)))
  
  # count total UMIs per cell
  umi_per_cell <- dt[, .N, by = cell]
  
  # visualize umi per cell to determine lower bound of cells to filter
  ggplot(, aes(x = N)) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    labs(title = "UMIs per cell", x = "number of UMIs (log10)", y = "cell count")
  
  # filter out low UMI cells (less than 500)
  # (could point to empty droplets or dead/broken cells)
  filtered_cells <- umi_per_cell[N >= 500, cell] # this outputs a vector of the barcode indices
  dt_filtered <- dt[cell %in% filtered_cells] # grabs those indices and selects from datatable
  
  # now get the barcodes that match the filtered indices
  filtered_barcodes <- barcodes[filtered_cells]
  
  # reset the cell index based on the filtering by adding a new column
  cell_index_map <- setNames(seq_along(filtered_cells), filtered_cells) 
  dt_filtered[, cell_new := cell_index_map[as.character(cell)]]
  
  # remove any possible duplicates
  dt_unique <- unique(dt_filtered)
  
  # build the sparse UMI count matrix
  sparse_counts <- sparseMatrix(
    i = dt_unique$gene, # row indices
    j = dt_unique$cell_new, # column indices
    x = rep(1L, nrow(dt_unique)), # entries
    dims = c(length(gene_ids), length(filtered_barcodes)), # matrix dimensions
    dimnames = list(genes, filtered_barcodes)
  )
  
  # write files (matrix market directory)
  # seurat will reference these 3 files when constructing 
  # the counts matrix with Read10X()
  basename <- tools::file_path_sans_ext(file)
  directory <- paste0('filtered_feature_bc_matrix_', basename)
  matrix_file <- paste0(directory, "/matrix.mtx") # in matrix market format
  barcodes_file <- paste0(directory, "/barcodes.tsv") # list of barcodes (corresponds to matrix colnames)
  features_file <- paste0(directory, "/features.tsv") # list of gene ids (corresponds to matrix rownames)
  
  dir.create(directory, showWarnings = FALSE)
  writeMM(sparse_counts, file = matrix_file)
  writeLines(filtered_barcodes, con = barcodes_file)
  
  features_df <- data.frame(genes,
                            gene_ids,
                            feature_type)
  write.table(features_df,
              file = features_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # zip the files
  gzfile(matrix_file)
  gzfile(barcodes_file)
  gzfile(features_file)
}

# clear environment to conserve memory
rm(list=ls())