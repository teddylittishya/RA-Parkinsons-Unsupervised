
#' Generate a null data set by permutation
#'
#' Generates a null data set by either permuting the data set (therefore maintaining marginal
#' distributions but breaking associations across variables) or by sampling uniformly over a
#' hypercube (equivalent probabilities for categorical variables).
#'
#' Features detected to be all 0/1 are treated as binary categorical variables.
#'
#' @param x data set (data frame)
#' @import dplyr
#' @return 
sample_data_set = function(x, type=c("uniform", 'permuted')) {
    type = match.arg(type)
    if(type == 'uniform') {
        n = nrow(x)
        x = x %>% mutate(across(everything(), 
            function(.x) {
                .x = .x[!is.na(.x)]
                if(all(.x == 0 | .x == 1)) { # binary
                    return(sample(c(0, 1), size=n, replace=T))
                } else if(is.numeric(.x)) { # continuous
                    return(runif(n, min=min(.x), max=max(.x)))
                } else if(is.factor(.x)) { # categorical
                    return(factor(sample(unique(.x), size=n, replace=T)))
                } else { # other?
                    return(sample(.x, size=n, replace=T))
                }
            })) 
        return(x)
    } else if(type == 'permuted') {
        return(x %>% mutate(across(everything(), ~ sample(.x, size=nrow(x)))))
    }
}

#' Generate a random data set by random sampling from uniform distributions
#'
#' @param m number of rows
#' @param d number of columns
#' @param type one of `cont` or `cat` for continuous or categorical random data
#' @param num.cat number of categories for categorical random data; if length of num.cat is equal
#'        to d, then the number of categories for column i be num.cat[i] 
random_data_set = function(m, d, type=c('cont', 'cat'), num.cat=3) {
    type = match.arg(type)
    tmp = NULL
    if(type == 'cont') {
        tmp = do.call('data.frame', 
                      lapply(1:d, function(.x) runif(n=m)))
        colnames(tmp) = paste0('rand_cont_runif_', 1:n)
    } else if(type == 'cat') {
        if(length(num.cat) == d) {
            num.cat = num.cat
        } else {
            num.cat = rep(num.cat, d)
        }
        tmp = do.call('data.frame', 
                      lapply(num.cat, function(.x) factor(sample(1:.x, size=m, replace=T))))
        colnames(tmp) = paste0('rand_cat', num.cat, '_', 1:n)
    }
    return(tmp)
}



#' Hopkins statistic
#'
#' Calculate the random forest Hopkins statistic on a data set
#'
#' @param x data set
#' @param k kth nearest neighbor
#' @param m sample size
#' @param seed random seed
#' @param ntree number of trees
#' @return Hopkins statistic
Hopkins.rf = function(x, k=1, m=0.1*nrow(x), seed=42, ntree=1000, type=c("uniform", 'permuted')) {
    type = match.arg(type)
    if(!is.null(seed)) set.seed(seed)
    x_true = x %>% slice_sample(n=m)
    d = ncol(x_true)
    rf = randomForest(x=x_true, 
                      proximity=T, 
                      keep.forest=F, 
                      ntree=ntree)
    prox = rf$proximity
    diag(prox) = -1
    x_i = (1 - apply(prox, 1, function(.x) 
                     sort(.x, partial=(m - k + 1))[(m - k + 1)]))^(1.5 * log(d) * log(m))

    x_fake = sample_data_set(x, type=type) %>% slice_sample(n=m) 
    rf = randomForest(x=x_fake, 
                      proximity=T, 
                      keep.forest=F, 
                      ntree=ntree)
    prox = rf$proximity
    diag(prox) = -1
    y_i = (1 - apply(prox, 1, function(.x) 
                     sort(.x, partial=(m - k + 1))[(m - k + 1)]))^(1.5 * log(d) * log(m))

    sum(y_i) / (sum(x_i) + sum(y_i))
}


#' Hopkins statistic p-value
#'
#' Return a p-value for a given Hopkins statistic
#'
#' @param H Hopkins statistic (produced by `Hopkins.rf`)
#' @param m sample size used to calculated Hopkins statistic
#' @return numeric p-value
Hopkins.rf.pval = function(H, m) {
    pbeta(H, m, m, lower.tail=F)
}

#' Pairwise feature distance
#'
#' @param rf random forest
feature_dist = function(rf, mc.cores=1) {
    forest = rf$forest
    stopifnot(!is.null(forest))
    d = nrow(importance(rf)) # number of variables
    adj = Reduce('+', mclapply(1:rf$ntree, function(k) {
        ar = matrix(0, nrow=d, ncol=d)
        bestvar = forest$bestvar[,k]
        ld = forest$treemap[,1,k]
        rd = forest$treemap[,2,k]
        for(i in 1L:length(bestvar)) {
            if(bestvar[i] == 0) next
            if(bestvar[ld[i]] > 0) 
                ar[bestvar[i], bestvar[ld[i]]] = ar[bestvar[i], bestvar[ld[i]]] + 1
            if(bestvar[rd[i]] > 0) 
                ar[bestvar[i], bestvar[rd[i]]] = ar[bestvar[i], bestvar[rd[i]]] + 1
        }
        ar
    }, mc.cores=mc.cores))
    colnames(adj) = rownames(adj) = rownames(importance(rf))
    adj = adj + t(adj)

    fact = table(rf$forest$bestvar)
    fact = as.numeric(fact[as.character(1:nrow(adj))])
    fact[is.na(fact)] = 1
    t(adj / fact) / fact
}
