## This script performs "ICA", "IPCA", "NMF", "Diffusion Map", "MDS" and "AutoEncoder" for dimension reductions (DR) method comparison on a given scRNA dataset
## Next it incorporates the "PCA", "UMAP" & "tSNE" results from the PCA_UMAP_tSNE.R script and perform the Hierarchical & 
## K-means clustering for all 9 DR methods and store ARI & NMI scores along with the plots
# Author : Malki Wijesinghe
# Date : 2023-12-07

library(dplyr)
library(destiny)
library(keras)
library(tensorflow)
library(ggplot2)
library(aricode)
library(FactoMineR)
library(Seurat)
library(factoextra)
library(fpc)
library(SingleR)
library(gridExtra)
library(mixOmics)
library(formattable)
library(ggplot2)
library(ica)
library(RcppML)
library(kohonen)
library(psych)
library(Rtsne)


# This function is used to preprocess the data
preprocess_df <- function(df) {
  # Convert the sparse matrix to a regular matrix
  dense_matrix <- as(df, "matrix")
  # Convert to a data frame
  df_dense <- as.data.frame(dense_matrix, stringsAsFactors = FALSE)
  df <- df_dense %>% mutate_all(~ifelse(is.na(.), 0, .))
  df = df[rowSums(df != 0) > 0, ]
  # Normalize the data
  normalized_data <- NormalizeData(df) 
  df_normalized <- as.data.frame(t(normalized_data))
  return(df_normalized)
}


# This function is used to perform diffusion map
get_dm <- function(df_normalized) {
  set.seed(2023)
  res_dm <- DiffusionMap(df_normalized, n_eigs=5)
  res_dm <- eigenvectors(res_dm)[, seq_len(5), drop = FALSE]
  return (res_dm)
}


# This function is used to perform Independent principal component analysis
get_ipca <- function(df_normalized) {
  set.seed(2023)
  trans.ipca <- ipca(df_normalized, ncomp = 5) 
  return (trans.ipca$x)
}


# This function is used to perform Independent component analysis
get_ica <- function(df_normalized) {
  set.seed(2023)
  res_ica <- icafast(df_normalized, nc =5, maxit = 100, tol = 1e-6, alg = "par", fun = "exp", alpha = 1) 
  return (res_ica$S)
}


# This function is used to perform Non-negative matrix factorization
get_nmf <- function(df_normalized) {
  set.seed(2023)
  nmf_model <- RcppML::nmf(df_normalized, 100, tol = 1e-3)   
  return (nmf_model$w)
}


# This function is used to perform Multidimensional Scaling
get_mds <- function(df_normalized) {
  set.seed(2023)
  res_mds <- cmdscale(dist(df_normalized), k=5)
  return (res_mds)
}


# This function is used to perform Auto-encoder
get_ae <- function(df_normalized) {
  set.seed(2023)
  suppressPackageStartupMessages(library(keras))
  x_train = df_normalized
  # set training data
  x_train <- as.matrix(x_train)
  # set model
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 256, activation = "tanh", input_shape = ncol(x_train)) %>%
    layer_dense(units = 100, activation = "tanh", name = "bottleneck") %>%
    layer_dense(units = 256, activation = "tanh") %>%
    layer_dense(units = ncol(x_train))
  # view model layers
  model_summary <- summary(model)
  # compile model
  model %>% compile(loss = "mean_squared_error", optimizer = "adam")
  # fit model
  model %>% fit(x = x_train, y = x_train, epochs = 500, verbose = 0)
  # evaluate the performance of the model
  mse.ae2 <- print(paste0("MSE : ", evaluate(model, x_train, x_train)))
  # output <- c(model_summary, mse.ae2)
  # writeLines(output, con = file.path("C:/Users/malki/OneDrive/SFU Materials/Course Work/Fall 2023/CPSC 545/Lectures/cpsc545-project/results/AE"
  #                                    , file, "AE_summary.txt"))
  # extract the bottleneck layer
  intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
  intermediate_output <- predict(intermediate_layer_model, x_train)
  return (intermediate_output)
}



#Use below function to perform Hierarchical & K-means clustering on dimension reduced data and obtain ARI & NMI scores
get_clustering <- function(dr_output, base_data_path, file, method, tenx_labels, singleR_labels, output_path){
  #============ Clustering evaluation using 10X-Genomics graph based clusters as ground truth ============#
  labels = tenx_labels
  labels = as.numeric(labels$Cluster)
  
  # Perform hierarchical clustering
  k = as.numeric(length(unique(labels)))
  hc.res <- eclust(dr_output, "hclust", k = k , hc_metric = "euclidean", hc_method = "ward.D2", graph = FALSE)
  # Calculate Adjusted Rand Index (ARI) metrics
  clust_stats <- cluster.stats(d = dist(df_normalized), labels, hc.res$cluster)
  hc_ari1 <- clust_stats$corrected.rand
  # Calculate Normalized Mutual Information (NMI) for hierarchical clustering
  hc_nmi1 <- NMI(hc.res$cluster,labels)
  
  # Perform K-means clustering
  km.res <- eclust(dr_output, "kmeans", k = k, nstart = 25, graph = FALSE)
  # Calculate Adjusted Rand Index (ARI) metrics
  clust_stats <- cluster.stats(d = dist(df_normalized), labels, km.res$cluster)
  km_ari1 <- clust_stats$corrected.rand
  # Calculate Normalized Mutual Information (NMI) for K-means clustering
  km_nmi1 <- NMI(km.res$cluster,labels)
  
  
  #============ Clustering evaluation using singleR cell annotation as ground truth ============# 
  # Get the cell annotation using singleR
  labels = as.numeric(as.factor(singleR_labels))
  
  # Perform hierarchical clustering
  k = as.numeric(length(unique(labels)))
  hc.res <- eclust(dr_output, "hclust", k = k , hc_metric = "euclidean", hc_method = "ward.D2", graph = FALSE)
  # Calculate Adjusted Rand Index (ARI) metrics
  clust_stats <- cluster.stats(d = dist(df_normalized), labels, hc.res$cluster)
  hc_ari2 <- clust_stats$corrected.rand
  # Calculate Normalized Mutual Information (NMI) for hierarchical clustering
  hc_nmi2 <- NMI(hc.res$cluster,labels)
  
  # Perform K-means clustering
  km.res <- eclust(dr_output, "kmeans", k = k, nstart = 25, graph = FALSE)
  # Calculate Adjusted Rand Index (ARI) metrics
  clust_stats <- cluster.stats(d = dist(df_normalized), labels, km.res$cluster)
  km_ari2 <- clust_stats$corrected.rand
  # Calculate Normalized Mutual Information (NMI) for K-means clustering
  km_nmi2 <- NMI(km.res$cluster,labels)
  
  #============ Store ARI & NMI values in a character vector ============# 
  fileConn <- file(file.path(output_path, file, paste0(method,".txt")), "w")
  # Write ARI and NMI values to the file
  cat(paste("10X-Genomics ARI for hierarchical Clustering for ", hc_ari1), "\n", file = fileConn)
  cat(paste("10X-Genomics NMI for hierarchical Clustering for", hc_nmi1), "\n", file = fileConn)
  cat(paste("10X-Genomics ARI for K-means Clustering for ", km_ari1), "\n", file = fileConn)
  cat(paste("10X-Genomics NMI for K-means Clustering for",  km_nmi1), "\n", file = fileConn)
  cat(paste("singleR ARI for hierarchical Clustering for ", hc_ari2), "\n", file = fileConn)
  cat(paste("singleR NMI for hierarchical Clustering for",  hc_nmi2), "\n", file = fileConn)
  cat(paste("singleR ARI for K-means Clustering for ", km_ari2), "\n", file = fileConn)
  cat(paste("singleR NMI for K-means Clustering for", km_nmi2), "\n", file = fileConn)
  # Close the file connection
  close(fileConn)
}



# This function generates the tSNE plots after each dimension reduction method
get_tsne <- function(dr_output, method, base_data_path, tenx_labels, singleR_labels) {
  tsne_output <- Rtsne(dr_output, dims = 2, perplexity = 30, verbose = TRUE, check_duplicates = FALSE)
  # Convert t-SNE output to a data frame
  tsne_df <- as.data.frame(tsne_output$Y)
  colnames(tsne_df) <- c("t_SNE_1", "t_SNE_2")
  
  #Get tSNE plot with 10xgenomics cluster labels
  labels_10XGenomics = as.factor(tenx_labels)
  
  tsne_plot1 <- ggplot(tsne_df, aes(x = t_SNE_1, y = t_SNE_2, color = labels_10XGenomics)) +
    geom_point(size = 1) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = method) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 20, hjust = 0.5)  # Center the plot title
    )
  
  #Get tSNE plot with singleR annotated cell labels
  tsne_plot2 <- ggplot(tsne_df, aes(x = t_SNE_1, y = t_SNE_2, color = singleR_labels)) +
    geom_point(size = 1) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = method) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 20, hjust = 0.5)  # Center the plot title
    )
  return(list(tsne_plot1,tsne_plot2))
}


# This function is used to map each dimension reduction output to the relevant method names
assign_dr_output <- function(method) {
  if (method == "ICA") {
    dr_output = dr_output_ica
  } else if (method == "IPCA") {
    dr_output = dr_output_ipca
  } else if (method == "Diffusion Map") {
    dr_output = dr_output_dm
  } else if (method == "NMF") {
    dr_output = dr_output_nmf
  } else if (method == "MDS") {
    dr_output = dr_output_mds
  } else if (method == "AutoEncoder") {
    dr_output = dr_output_ae
  } else if (method == "tSNE") {
    dr_output = dr_output_tsne
  } else if (method == "UMAP") {
    dr_output = dr_output_umap
  } else if (method == "PCA") {
    dr_output = dr_output_pca
  } else {
    print("Error - Check dimension reduction method names")
  }
  return(dr_output)
}





########################################################################## Execution ##########################################################################

base_data_path = "C:/Users/malki/OneDrive/SFU Materials/Course Work/Fall 2023/CPSC 545/Lectures/cpsc545-project/data"
output_path = "C:/Users/malki/OneDrive/SFU Materials/Course Work/Fall 2023/CPSC 545/Lectures/cpsc545-project/results"
files = list("Brain_Tumor_3p_LT_filtered_feature_bc_matrix", "Breast_Cancer_3p_LT_filtered_feature_bc_matrix", "Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer_filtered_feature_bc_matrix")
methods = list("tSNE","UMAP","PCA","ICA", "IPCA", "NMF", "Diffusion Map", "MDS","AutoEncoder")


for (file in files) {
  #=============================== Data Pre-processing ===============================#
  #  Read the scRNA data set
  df <- Read10X_h5(filename = file.path(base_data_path, paste0(file,".h5")))
  print(dim(df))
  df_normalized = preprocess_df(df)
  
  # Get the singleR cell annotation
  ref <- celldex::BlueprintEncodeData()
  pred <- SingleR(test = df_normalized, ref = ref, labels = ref$label.main)
  singleR_labels = pred$labels
  # Import the 10xGenomics cluster labels
  labels = read.csv(file.path(base_data_path, paste0(file,".csv")))
  tenx_labels = as.factor(labels$Cluster)
  
  # Get the t-SNE plots for Pre-processed data for visualization with 10xGenomics labels & singleR cell annotation
  tsne_output <- Rtsne(as.matrix(df_normalized), dims = 2, perplexity = 30, verbose = TRUE)
  tsne_df <- as.data.frame(tsne_output$Y)
  colnames(tsne_df) <- c("t_SNE_1", "t_SNE_2")
  ## 10xGenomics labels
  tsne_plot_df1 <- ggplot(tsne_df, aes(x = t_SNE_1, y = t_SNE_2, color = tenx_labels)) +
    geom_point(size = 1) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = "Preprocessed Dataset") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 20, hjust = 0.5),  # Center the plot title
      legend.position ="right", legend.text = element_text(size = 12)
    )
  ## singleR annotation
  tsne_plot_df2 <- ggplot(tsne_df, aes(x = t_SNE_1, y = t_SNE_2, color = singleR_labels)) +
    geom_point(size = 1) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = "Preprocessed Dataset") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 20, hjust = 0.5),  # Center the plot title
      legend.position ="right", legend.text = element_text(size = 12)
    )
  
  # Extract the legends to be used in the grid plots
  legend1 <- cowplot::get_legend(tsne_plot_df1 + theme(legend.position = "right"))
  legend2 <- cowplot::get_legend(tsne_plot_df2 + theme(legend.position = "right"))
  tsne_plot_df1 = tsne_plot_df1 + theme(legend.position="none")
  tsne_plot_df2 = tsne_plot_df2 + theme(legend.position="none")
  
  
  #=============================== Dimension Reduction & Evaluation ===============================#
  
  # Import the data matrices produced by tSNE, UMAP & PCA dimension reduction methods for each scRNA dataset
  dr_output_tsne = readRDS(file.path(output_path, "Riya", file, "tsne.rds" ))
  dr_output_umap = readRDS(file.path(output_path, "Riya", file, "umap.rds" ))
  dr_output_pca = readRDS(file.path(output_path, "Riya", file, "pca.rds" ))
  
  # Perform Dimension Reduction for the ICA, IPCA, NMF, Diffusion Map, MDS, AutoEncoder methods 
  dr_output_dm = get_dm(df_normalized)
  dr_output_ipca = get_ipca(df_normalized)
  dr_output_ica = get_ica(df_normalized)
  dr_output_nmf = get_nmf(df_normalized)
  dr_output_mds = get_mds(df_normalized)
  dr_output_ae = get_ae(df_normalized)
  
  # Create lists designated to hold the tSNE plots for all 9 dimension reduction techniques
  tsne_plots1 <- list() # for the tSNE plots with 10xGenomics labels
  tsne_plots1[["tsne_plot_df1"]] <- tsne_plot_df1
  tsne_plots2 <- list() # for the tSNE plots with singleR annotation
  tsne_plots2[["tsne_plot_df2"]] <- tsne_plot_df2
  
  # Iterate over each method, generating clustering results and creating visual representations for each.
  for (method in methods) {
    dr_output = assign_dr_output(method)
    # Clustering evaluation
    get_clustering(dr_output, base_data_path, file, method, tenx_labels, singleR_labels, output_path)
    # Visualization
    tsne_plot1_tsne_plot2 = get_tsne(dr_output, method, base_data_path, tenx_labels, singleR_labels)
    tsne_plot1 <- tsne_plot1_tsne_plot2[[1]]
    tsne_plot2 <- tsne_plot1_tsne_plot2[[2]]
    tsne_plot1 = tsne_plot1 + theme(legend.position="none")
    tsne_plot2 = tsne_plot2 + theme(legend.position="none")
    tsne_plots1[[method]] <- tsne_plot1
    tsne_plots2[[method]] <- tsne_plot2
  }
  
  #Save plots
  pdf(file.path(output_path, "aggregated_plots", file, "10x before & after DR.pdf"), width = 25, height = 10)
  cowplot::plot_grid(grid.arrange(grobs = tsne_plots1, ncol = 5), legend1, rel_widths = c(1, 0.08))
  dev.off()
  
  pdf(file.path(output_path, "aggregated_plots", file, "singleR before & after DR.pdf"), width = 25, height = 10)
  cowplot::plot_grid(grid.arrange(grobs = tsne_plots2, ncol = 5), legend2, rel_widths = c(1, 0.08))
  dev.off()
}
