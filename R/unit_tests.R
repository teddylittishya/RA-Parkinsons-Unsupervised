source("cluster_module.R")
source("evaluation_module.R")

#' Test the cluster module and evaluation module on the iris data set (data frame)
#' @param ... params to pass to `cluster_module`
iris.test1 = function(...) {
    cvec = cluster_module(iris[,1:4], ...)
    eval = evaluation_module(cvec, "hi")
    list(cluster_sizes=table(cvec), evaluation_metrics=eval)
}

#' Test the cluster module and evaluation module on the iris data set (dist)
#' @param ... params to pass to `cluster_module`
iris.test2 = function(method="euclidean", ...) {
    iris.dist = dist(iris[, 1:4], method=method)
    cvec = cluster_module(iris.dist, ...)
    eval = evaluation_module(iris.dist, cvec, "hi")
    list(cluster_sizes=table(cvec), evaluation_metrics=eval)
}

result2 = iris.test2(method="euclidean",
                     algo = "kmeans",
                     k = 10)
result2$cluster_sizes
result2$evaluation_metrics

result2b = iris.test2(method="euclidean",
                     algo = "pam",
                     k = 10)
result2b$cluster_sizes
result2b$evaluation_metrics

result2c = iris.test2(method="euclidean",
                     algo = "fanny",
                     k = 10)
result2c$cluster_sizes
result2c$evaluation_metrics

result3 = iris.test2(method="euclidean",
                     algo = "hierarchical",
                     agg_method = "average",
                     k = 10)
result3$cluster_sizes
result3$evaluation_metrics

result3 = iris.test2(method="euclidean",
                     algo = "agnes",
                     agg_method = "average",
                     k = 10)
result3$cluster_sizes
result3$evaluation_metrics

result3 = iris.test2(method="euclidean",
                     algo = "diana",
                     k = 10)
result3$cluster_sizes
result3$evaluation_metrics

result4 = iris.test2(method="euclidean",
                     algo = "OPTICS",
                     minPts = 10,
                     xi = 0.001)
result4$cluster_sizes
result4$evaluation_metrics

