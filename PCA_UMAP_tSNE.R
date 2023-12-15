# This is the script to perform PCA, UMAP and tSNE for the datasets using Seurat R package
# The final output expected is QC plot, normalized counts matrix and dimensionally reduced matrices which will be used for clustering in dim_reduction_and_clustering.R script
# Author : Riya Saju
# Date : 2023-12-07

#load libraries
library(Seurat)
library(tidyverse)
library(hdf5r)
library(SingleR)
library(celldex)
library(pheatmap)
library(gridGraphics)

#function to process datasets
process_dataset <- function(dataset_filename, output_prefix) {
  
  # 1. Load the dataset
  dataset <- Read10X_h5(filename = dataset_filename)
  str(dataset)
  # Convert the sparse matrix to a regular matrix 
  dense_matrix <- as(dataset, "matrix")
  # Convert to a data frame
  df_dense<- as.data.frame(dense_matrix, stringsAsFactors = FALSE)
  df <- df_dense %>% mutate_all(~ifelse(is.na(.), 0, .))
  df = df[rowSums(df != 0) > 0, ]
  # Display the original data
  print("Original Data:")
  print(head(df))
  
  # 2. Initialize the Seurat object with the raw (non-normalized data)
  seurat.obj <- CreateSeuratObject(counts = df)
  #View(seurat.obj@meta.data)
  
  # 3. QC of datasets
  #Not using this for now - no much difference, can't reduce cells further(very less no. of cells),already using filtered dataset
  #seurat.obj$mitoPercent <-PercentageFeatureSet(seurat.obj, pattern = '^MT-')
  #seurat.obj.filtered <- subset(seurat.obj, subset = mitoPercent < 10)
  FeatureScatter(seurat.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')+geom_smooth(method = 'lm')+theme(text = element_text(size = 8))
  ggsave(filename = paste0("/Users/riyasaju/Downloads/CPSC545/", output_prefix, "_qcplot.pdf"), plot = last_plot())
  
  # 4. Normalize data
  seurat.obj <- NormalizeData(seurat.obj)
  #Creating a df for clustering 
  df_normalized <- NormalizeData(df)
  # Convert to a data frame and saving as an RDS for clustering 
  matrix <- as.data.frame(t(df_normalized))
  saveRDS(matrix, file = paste0("/Users/riyasaju/Downloads/CPSC545/",output_prefix, "_normalizedcounts.rds"))
  
  # 5. Identify highly variable features
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  
  # 6. Scaling
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  str(seurat.obj)
  
  # 7. Perform Linear dimensionality reduction
  #Perform PCA
  seurat.obj<- RunPCA(seurat.obj)
  #Extracting the PCA results
  pca_results<-seurat.obj@reductions[["pca"]]
  # View the PCA results
  head(pca_results)
  # Extract the PCA matrix of PCA results and save as an RDS for clustering script
  pca_matrix<-pca_results@cell.embeddings
  saveRDS(pca_matrix, file = paste0("/Users/riyasaju/Downloads/CPSC545/",output_prefix, "_pca.rds"))
  
  # 8. Perform non-linear dimensionality reduction 
  #Perform umap
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)
  #Extracting the UMAP results
  umap_results<-seurat.obj@reductions[["umap"]]
  # Extract the UMAP matrix of UMAP results and save as an RDS for clustering script
  umap_matrix<-umap_results@cell.embeddings
  saveRDS(umap_matrix, file = paste0("/Users/riyasaju/Downloads/CPSC545/",output_prefix, "_umap.rds"))
  
  # Perform t-SNE
  seurat.obj <- RunTSNE(seurat.obj, dims = 1:15)
  #Extracting the TSNE results
  tsne_results<-seurat.obj@reductions[["tsne"]]
  # Extract the TSNE matrix of TSNE results and save as an RDS for clustering script
  tsne_matrix<-tsne_results@cell.embeddings
  saveRDS(tsne_matrix, file = paste0("/Users/riyasaju/Downloads/CPSC545/",output_prefix, "_tsne.rds"))
}

#function for looping the process_dataset function on a directory of input datasets
process_datasets_in_folder <- function(folder_path) {
  # Get a list of file names in the folder
  file_names <- list.files(path = folder_path, pattern = "\\.h5$", full.names = TRUE)
  
  # Loop over each file and process the dataset
  for (file_name in file_names) {
    # Extract the dataset name from the file name
    dataset_name <- tools::file_path_sans_ext(basename(file_name))
    
    # Call the processing function
    process_dataset(file_name, dataset_name)
  }
}

# Running the analysis using the functions created
process_datasets_in_folder("/Users/riyasaju/Downloads/CPSC545/datasets/")

