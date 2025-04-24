#' Clustering module
#'
#' Takes as input a distance matrix and runs clustering algorithms with a given set of parameters.
#'
#' @param x distance matrix object of type `dist` with dimensions n x n
#' @param algo clustering algorithm
#' @param ... additional arguments depending on the choice of `algo`
#'  - For 'OPTICS': 
#'     \itemize{
#'       \item \code{minPts}: Minimum number of points in a neighborhood for a point to be considered a core point.
#'       \item \code{xi}: Minimum reachability distance for clustering.
#'     }
#'   - For 'kmeans': 
#'     \itemize{
#'       \item \code{k}: Number of clusters or initial cluster centers.
#'     }
#'   - For 'hierarchical': 
#'     \itemize{
#'       \item \code{agg_method}: Linkage method ('complete', 'single', 'average', etc.) for hierarchical clustering.
#'       \item \code{k}: Number of clusters to be formed by cutting the dendrogram.
#'     }
#' @return vector of length n with cluster memberships
#' @import dbscan
#' @import cluster
#' @import mclust
#' @export
#' 
cluster_module <- function(x, 
                           algo=c('OPTICS', 'kmeans', 'pam', 'fanny', 'hierarchical', 'agnes', 'diana','spectral', 'affinity', 'GMM'),
                           ...) {
  algo = match.arg(algo)
  params <- list(...)
  cvec <- vector()
  switch(algo,
         'OPTICS' = {
           # Ordering Points to Identify the Clustering Structure (OPTICS)
           reach <- dbscan::optics(x, minPts=params$minPts)
           clusters <- dbscan::extractXi(reach, xi=params$xi, minimum=T)
           cvec <- clusters$cluster
         },
         'kmeans' = {
           # K-Means Clustering (Partitioning Clustering)
           clusters <- stats::kmeans(x, centers = params$k)
           cvec <- clusters$cluster
         },
         'pam' = {
           # Partitioning Around Medoids (Partitioning Clustering)
           clusters <- cluster::pam(x, k = params$k)
           cvec <- clusters$clustering
         },
         'fanny' = {
           # Fuzzy Analysis Clustering (Partitioning Clustering)
           clusters <- cluster::fanny(x, k = params$k, maxit=1500)
           cvec <- clusters$clustering
         },
         'hierarchical' = {
           # Hierarchical Clustering
           clusters <- stats::hclust(x, method = params$agg_method)
           cvec <- stats::cutree(clusters, k = params$k)
         },
         'agnes' = {
           # Agglomerative Nesting (Hierarchical Clustering)
           if(params$agg_method == 'mcquitty') params$agg_method = 'weighted'
           clusters <- cluster::agnes(x, method=params$agg_method)
           cvec <- stats::cutree(clusters, k = params$k)
         },
         'diana' = {
           # Divisive ANAlysis Clustering (Hierarchical Clustering)
           clusters <- cluster::diana(x)
           cvec <- stats::cutree(clusters, k = params$k)
         },
         'spectral' = {
           # Spectral Clustering
           clusters <- kernlab::specc(as.matrix(x), centers = params$k)
           cvec <- as.numeric(clusters)
         },
         'affinity' = {
           # Affinity Propagation
           clusters <- apcluster::apcluster(negDistMat(r=2), as.matrix(x))
           cvec <- clusters@idx
         },
         'GMM' = {
           # Gaussian Mixture Model Clustering (Mclust)
           clusters <- mclust::Mclust(as.matrix(x), G = params$G)
           cvec <- clusters$classification
         },
         {
           stop("The argument `algo` is required. ", 
                "Possible values are ",
                'OPTICS', 'kmeans', 'pam', 'fanny', 'hierarchical', 'agnes', 'diana','spectral', 
                'affinity', 'GMM')
         })
  return(cvec)
}



