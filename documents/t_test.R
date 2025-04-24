

# Source the main clustering module
source("./R/cluster_module.R")

# Load necessary libraries
library(randomForest)
library(cluster)
library(dbscan)
library(stats)
library(mclust)
library(kernlab)
library(apcluster)

# Prepare the data and distance matrix
data(iris)
iris_data <- iris[, -5]  # Remove the species column for clustering
rf_model <- randomForest(iris_data, ntree = 500, proximity = TRUE)
distance_matrix <- as.dist(1 - rf_model$proximity)

# Test each algorithm
cat("Testing different clustering algorithms...\n")

# OPTICS
cat("\nOPTICS:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'OPTICS', minPts = 5, xi = 0.05)
print(cluster_result)
print(table(cluster_result, iris$Species))

# K-Means
cat("\nK-Means:\n")
cluster_result <- cluster_module(iris_data, algo = 'kmeans', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# PAM (Partitioning Around Medoids)
cat("\nPAM:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'pam', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# FANNY (Fuzzy Analysis Clustering)
cat("\nFANNY:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'fanny', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# Hierarchical Clustering
cat("\nHierarchical Clustering:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'hierarchical', agg_method = 'average', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# AGNES (Agglomerative Nesting)
cat("\nAGNES:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'agnes', agg_method = 'average', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# DIANA (Divisive Analysis)
cat("\nDIANA:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'diana', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# Spectral Clustering
cat("\nSpectral Clustering:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'spectral', k = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

# Affinity Propagation
cat("\nAffinity Propagation:\n")
cluster_result <- cluster_module(distance_matrix, algo = 'affinity')
print(cluster_result)
print(table(cluster_result, iris$Species))

# Gaussian Mixture Model (GMM)
cat("\nGaussian Mixture Model (GMM):\n")
cluster_result <- cluster_module(iris_data, algo = 'GMM', G = 3)
print(cluster_result)
print(table(cluster_result, iris$Species))

