#library(fpc)
#library(clusterSim)
# These import calls don't belong in packages. import them separately and within functions 
# use the package::function syntake


#' Evaluation module
#'
#' Takes as input the cluster memberships produced by a single run of `cluster_module`, and the
#' parameters used by `cluster_module`. Returns a list of evaluation metrics for that run of
#' `cluster_module`.
#'
#' @param x distance matrix
#' @param cvec numeric vector of length n with cluster memberships
#' @param parameter list of parameters used by `cluster_module`
#' @param supervised boolean indicating whether or not to use supervised learning as a metric of
#'   cluster quality
#' @return numeric vector of evaluation metrics including
#'   - Cohesiveness: Measures how tightly grouped points are within clusters.
#'   - Separation: Measures how distinct clusters are from each other.
#'   - Cluster Sizes: Number of data points in each cluster.
#'   - Other relevant clustering metrics depending on the implementation.
#' @import fpc
#' @import clusterSim
#' @export
#' 
evaluation_module <- function(x,
                              cvec, 
                              parameter, 
                              supervised=F) {
  # calculate clustering statistics
  stats <- fpc::cluster.stats(x, cvec)
  if(is.null(cvec) | all(cvec == 0)) { 
    return(list(nbclust=NA, cohesiveness=NA, distinctiveness=NA, n_between=NA, n_within=NA,
                avg_silwidth=NA, max_diameter=NA, min_separation=NA, pearson_gamma=NA,
                dunn_index=NA, dunn2_index=NA, entropy=NA, wb_ratio=NA, ch_index=NA, widest_gap=NA,
                sep_index=NA, index_C=NA, index_G2=NA, index_G3=NA, index_S=NA
    )) 
  }
  # if(algo == "OPTICS")
  which.zero = which(cvec == 0)
  cvec = cvec[!which.zero]
  x = x[!which.zero, !which.zero]
  cvec = as.numeric(factor(cvec)) # renumber cvec
  
  # prepare output with evaluation metrics
  out = list(nbclust = stats$cluster.number,
             cohesiveness = stats$average.within, 
             # Measures how closely points are grouped within the same cluster.
             distinctiveness = stats$average.between, 
             # Measures how well clusters are separated from each other.
             n_between = stats$n.between,
             # Indicates the number of between-cluster distance measurements.
             n_within = stats$n.within,
             # Indicates the number of within-cluster distance measurements.
             avg_silwidth = stats$avg.silwidth,
             #  Provides the average silhouette width, which measures the quality of the clusters.
             max_diameter = stats$max.diameter,
             # Measures the greatest distance within two points in a cluster
             min_separation = stats$min.separation,
             # Measures the minimum separation between clusters.
             # gamma_coef = stats$g2,
             # Represents the gamma coefficient, used to evaluate clustering quality.
             # g3_coeff = stats$g3,
             # Represents the g3 coefficient, used to evaluate clustering quality.
             pearson_gamma = stats$pearsongamma,
             # Represents the pearson gamma coefficient which assess clustering quality. 
             dunn_index = stats$dunn,
             # Evaluates the ratio of minimum inter-cluster distance to maximum intra-cluster distance.
             dunn2_index = stats$dunn2,
             # Measures the modified Dunn index.
             entropy = stats$entropy,
             # Measures the entropy (distance and homogenity) of the clusters
             wb_ratio = stats$wb.ratio,
             # Represents the within-between cluster ratio, used for clustering evaluation.
             ch_index = stats$ch,
             # Provides the Calinski-Harabasz index, used to evaluate the dispersion of clusters.
             widest_gap = stats$widestgap,
             # Measures the widest gap between clusters.
             sep_index = stats$sindex,
             # Represents the separation index, used for clustering evaluation.
             index_C = clusterSim::index.C(x, cvec),
             # Calculates Hubert & Levin C index - internal cluster quality index
             index_G2 = clusterSim::index.G2(x, cvec),
             # Calculates G2 internal cluster quality index
             index_G3 = clusterSim::index.G3(x, cvec),
             # Calculates G3 internal cluster quality index
             index_S = clusterSim::index.S(x, cvec)
             # Calculates Rousseeuw’s Silhouette internal cluster quality index
  )
  
  #assign class to output list
  class(out) = 'evaluation_module'
  return(out)
}
