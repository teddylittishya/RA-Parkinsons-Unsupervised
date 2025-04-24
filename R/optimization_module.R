#' Optimization module
#'
#' Takes as input a distance matrix and determines the best clustering algorithm and parameters.
#'
#' @param data distance matrix
#' @param parameters list of parameters
#' @param optimization_method optimization algorithm
#' @param evaluation_metric evaluation metric
#' @export
#'
optimization_module <- function(data
                                parameters,
                                optimization_method=c('grid', 'stgd', 'adam'), 
                                evaluation_metric=c()) {
    optimization_method = match.arg(optimization_method)
    evaluation_metric = match.arg(evaluation_metric)
    #optimg:: 
    out = list(
        optimization_method = optimization_method,
        evaluation_metric = evaluation_metric,
        optimization_results = data.frame(),
        best_params = list(),
        best_metrics = 1,
    )
    class(out) = 'optimization_module'
    return(out)
}
