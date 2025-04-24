# Implementation of the RCE algorithm

#' Consensus clustering for feature selection
#'
#' using `kmeans`
#' @param sD full data set
#' @param nbiter number of bagging iterations
#' @param n `nrow(sD)` number of samples
#' @param k `ncol(sD)` number of features
#' @param nclusters the number of clusters for `kmeans`
#' @param p_bag size of random feature subspace (must be less than k)
consensus <- function(sD, nbiter, n, nclusters, k, p_bag, vfun=message) {
    M <- NULL # consensus matrix
    KK <- NULL
    imp <- matrix(0, n, k - 1) # importance matrix
    isoob <- matrix(0, n, 1)

    for (i in 1:nbiter) {
        if(i %% 50 == 0 & !is.null(vfun)) {
            vfun(paste('bagging iteration', i, 'of', nbiter))
        }
        # Bootstrap sample
        instances <- sample(1:n, n, replace = TRUE)

        # Randomly select features
        features <- sample(1:(k - 1), p_bag, replace = FALSE)

        # Project data and oob onto feature subset
        bag <- sD[instances, features, drop = FALSE]
        instances_unique <- unique(instances)
        instances_oob <- setdiff(1:n, instances_unique)
        oob <- sD[instances_oob, features, drop = FALSE]
        n_oob <- nrow(oob)

        # Clustering of bootstrap samples
        L <- kmeans(bag, nclusters, nstart = 1)
        K1 <- L$cluster

        K <- rep(0, n)
        K[instances] <- K1

        M1 <- matrix(0, n, nclusters)
        M1[cbind(1:n, K)] <- 1
        M <- cbind(M, M1) # update the consensus matrix M with M1
        KK <- rbind(KK, t(K)) # record the cluster assignments in KK.

        # Calculates the Euclidean distance between oob instances and the 
        # cluster centers of the bootstrap sample
        dist_cluster <- matrix(0, n_oob, nclusters)

        for (l in 1:nclusters) {
            if (sum(K1 == l) == 1) {
                base_centres <- bag[K1 == l, , drop = FALSE]
            } else {
                base_centres <- colMeans(bag[K1 == l, , drop = FALSE])
            }
            dist_cluster[, l] <- apply(oob, 1, function(x) sqrt(sum((x - base_centres)^2)))
        }

        # Determine the initial cluster assignments for oob instances.
        K_before <- apply(dist_cluster, 1, which.min)

        for (f in 1:p_bag) {
            oob_f <- oob[, f]
            # Permute the values each feature in oob data.
            oob_f_permute <- sample(oob_f)
            oob_modifed_f <- oob
            oob_modifed_f[, f] <- oob_f_permute

            # Recalculate the Euclidean distance and determines the new cluster assignments 
            dist_cluster <- matrix(0, n_oob, nclusters)

            for (l in 1:nclusters) {
                if (sum(K1 == l) == 1) {
                    base_centres <- bag[K1 == l, , drop = FALSE]
                } else {
                    base_centres <- colMeans(bag[K1 == l, , drop = FALSE])
                }
                dist_cluster[, l] <- apply(oob_modifed_f, 1, function(x) sqrt(sum((x - base_centres)^2)))
            }

            K_after <- apply(dist_cluster, 1, which.min)

            # Increment the importance score if the cluster assignment changes.
            for (j in 1:n_oob) {
                if (K_after[j] != K_before[j]) {
                    variable <- features[f]
                    l <- instances_oob[j]
                    imp[l, variable] <- imp[l, variable] + 1
                }
            }
        }
    }

    return(list(M = M, imp = imp))
}

# Select features based on their importance scores
scree_test <- function(imp_var, ...) {
    if(length(dim(imp_var)) < 2) { imp_var <- t(imp_var) }
    k = ncol(imp_var)
    nclusters = nrow(imp_var)
    # Sort features by importance
    sorted_imp <- t(apply(imp_var, 1, sort, decreasing = TRUE))
    indices <- t(apply(imp_var, 1, order, decreasing = TRUE))

    df <- matrix(0, nclusters, k - 2) # store the differences between consecutive importance scores
    acc <- matrix(0, nclusters, k - 3) # store the differences of df differences
    sc <- matrix(0, nclusters, k - 4) # store the sum of the absolute values of consecutive acc 
                                      # differences.

    # Calculate differences for each cluster
    for (i in 1:nclusters) {
        for (j in 1:(k - 2)) {
            df[i, j] <- sorted_imp[i, j] - sorted_imp[i, j + 1]
        }
    }

    # Calculate acceleration differences for each cluster
    for (i in 1:nclusters) {
        for (j in 1:(k - 3)) {
            acc[i, j] <- df[i, j] - df[i, j + 1]
        }
    }

    # Calculate the scree values by finding the point where the sum of consecutive 
    # second-order differences (sc) is maximum
    scree <- rep(0, nclusters) 
    # stores the index of the feature at the maximum point for each cluster
    for (i in 1:nclusters) {
        max_val <- 0
        for (j in 1:(k - 4)) {
            sc[i, j] <- abs(acc[i, j]) + abs(acc[i, j + 1])
            if (sc[i, j] > max_val) {
                scree[i] <- j
                max_val <- sc[i, j]
            }
        }
    }

    # Select the top features up to the scree point for each cluster 
    selected_features <- matrix(0, nclusters, max(scree))
    sf <- c() # collects all selected feature indices across clusters

    for (i in 1:nclusters) {
        for (j in 1:scree[i]) {
            selected_features[i, j] <- indices[i, j]
            sf <- c(sf, indices[i, j])
        }
    }

    sf_unique <- unique(sf) # stores the unique indices of all selected features

    return(list(selected_features = selected_features, sf_unique = sf_unique))
}

#' RCE algorithm
#'
#' Takes as input a data matrix or dataframe and outputs the estimated local feature's importance
#' and the selected features.
#'
#' @param data data frame or matrix 
#' @param ncluster number of clusters
#' @param nbiter number of committee members
#' @param vfun verbose function
#' @return list of feature's importance and selected features
#' @export
#' 
run_RCE = function(data, nclusters, nbiter, vfun=message) {
    data = as.data.frame(data)
    k <- ncol(data) # Total number of features
    n <- nrow(data) # Number of instances
    p_bag <- round(sqrt(k)) # number of features for each bag

    # For normalization
    # Yi - I removed the caret dependency by manually range scaling the data
    data <- as.data.frame(apply(data, 2, function(x) { 
                        (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T)) 
    }))

    # Performs consensus clustering for feature selection
    consensus_result <- consensus(data, nbiter, n, nclusters, k + 1, p_bag, vfun=vfun)
    M <- consensus_result$M
    imp <- consensus_result$imp

    # Calculate the similarity matrix from the co-association matrix
    SK <- M %*% t(M)
    # Convert to dissimilarity matrix
    Diss <- 1 - SK / nbiter
    diag(Diss) <- 0
    # Perform hierarchical clustering
    Z <- hclust(as.dist(Diss), method = "average")
    K_before <- cutree(Z, k = nclusters)

    # Normalize feature importance scores
    for (i in 1:n) {
        for (j in 1:(k)) {
            imp[i, j] <- imp[i, j] / nbiter
        }
    }

    # Aggregate feature importance scores for each cluster
    imp_var <- matrix(0, nclusters, k)
    for (i in 1:n) {
        for (j in 1:(k)) {
            imp_var[K_before[i], j] <- imp_var[K_before[i], j] + imp[i, j]
        }
    }

    # Return the selected features
    selected_features_result <- scree_test(imp_var, nclusters, k)
    selected_features <- selected_features_result$selected_features
    sf_unique <- selected_features_result$sf_unique

    return(list(imp_var = imp_var, 
                selected_features = sf_unique,
                final_clustering = K_before,
                imp_samp = imp,
                imp = apply(imp, 2, mean)
                ))
}

#' RCE original
#'
#' Original function used by the RCE authors - takes as input a path to a data set and
#' automatically runs RCE. The true cluster labels are provided in the data set, which is not what
#' we want.
#'
#' @param data_path path of the data set
#' @param nbiter number of committee members
#' @param vfun verbose function
#' @return list of feature's importance and selected features
#' @export
#' 
RCE <- function(data_path, nbiter, vfun=message) {
  # Read the data
  data <- read.table(data_path, header = FALSE, sep = "\t")
  k <- ncol(data) # Total number of features
  n <- nrow(data) # Number of instances
  V <- data[, k] # True clustering
  nclusters <- max(V) # Number of clusters
  
  # Database without label
  data <- data[, 1:(k - 1)]

  run_RCE(data, nclusters = nclusters, nbiter = nbiter, vfun = vfun)
}

