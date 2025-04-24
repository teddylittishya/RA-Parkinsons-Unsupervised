run_all = function(data, 
                   labels,
                   k=seq(2, 10), 
                   agg = c('average','single','complete','ward','mcquitty'),
                   xi=c(0.01, 0.05, 0.1, 0.2, 0.5), 
                   minPts=4:10) {
  data.dist = dist(data, method="euclidean")
  #data.rf <- randomForest::randomForest(x=data, proximity=TRUE, keep.forest=FALSE, ntree=2000)
  #data.dist = 1 - data.rf$proximity
  
  
  algos = eval(formals(args(cluster_module))$algo)
  df1 = lapply(k, function(kval) {
    lapply(algos, function(algo) {
      if(algo == "OPTICS" | algo == "hierarchical" | algo == "agnes") { return(NULL) }
      cvec = cluster_module(data.dist, algo=algo, k=kval, agg_method="complete")
      names(cvec) = paste0("sample",1:length(cvec))
      c(param.algo=algo, param.k=kval, evaluation_module(data.dist, cvec), 
        NMI=aricode::NMI(unlist(cvec), labels),
        ARI=mclust::adjustedRandIndex(unlist(cvec), labels),
        cvec)
    })
  })
  df1 = bind_rows(unlist(df1, recursive=F))
  
  optics_params = expand.grid(xi=xi, minPts=minPts)
  df2 = apply(optics_params, 1, function(row) {
    algo = "OPTICS"
    xival = row[["xi"]]
    mpval = row[["minPts"]]
    cvec = cluster_module(data.dist, algo="OPTICS", xi=xival, minPts=mpval)
    if(is.null(cvec)) cvec = rep(0, nrow(data.dist))
    names(cvec) = paste0("sample",1:length(cvec))
    c(param.algo="OPTICS", param.xi=xival, param.minPts=mpval, 
      evaluation_module(data.dist, cvec), 
      NMI=aricode::NMI(unlist(cvec), labels),
      ARI=mclust::adjustedRandIndex(unlist(cvec), labels),
      cvec)
  })
  df2 = bind_rows(df2)
  
  
  df3 = bind_rows(lapply(k, function(kval) {
    bind_rows(lapply(algos, function(algo) {
      bind_rows(lapply(agg, function(aggval) {
        if (algo == "hierarchical" | algo == "agnes") {
          cvec = cluster_module(data.dist, algo=algo, k=kval, agg_method=aggval)
          names(cvec) = paste0("sample",1:length(cvec))
          c(param.algo=algo, param.k=kval, param.agg_method=aggval, 
            evaluation_module(data.dist, cvec), 
            NMI=aricode::NMI(unlist(cvec), labels),
            ARI=mclust::adjustedRandIndex(unlist(cvec), labels),
            cvec)
        }
      }))
    }))
  }))
  
  # eval <- bind_rows(df1, df2, df3) %>% relocate(starts_with('param'))
  eval <- bind_rows(df1, df3) %>% relocate(starts_with('param'))
  # cvec_mat <- do.call(rbind, cvec_val)
  # out = list(eval.res=eval, cvec_val.res=cvec_val)
  # class(out) = 'run_module'
  return(eval)
}