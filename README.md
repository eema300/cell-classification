# cell-classification

## About
The h5_analyze_demo.R script uses the [RegenOrNoRegen](https://github.com/neurohugo/RegenOrNoRegen) R package 
to classify neurons as regenerating or non-regenerating based on gene expression levels from scRNA-seq data.

## Usage
This is just a demo mostly for
* extracting from an h5 format that is incompatible with Seurat's `Read10X_h5()` function
* using the RegenOrNoRegen package
<br>
The script order is: h5_extract_demo.R -> h5_filter_options_demo.R -> h5_analyze_demo.R

## Dependencies
See the [RegenOrNoRegen](https://github.com/neurohugo/RegenOrNoRegen) repository for package dependency requirements.
