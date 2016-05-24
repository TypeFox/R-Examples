summary.traj <-
function(object, round.pos = 2, ...)
{
  cat("Number of observations: ", nrow(object$data),"\n\n")
  
  cat("Number of clusters: ", ncol(object$clust.distr),"\n\n")
  
  cat("Cluster distribution:\n")
  print(round(object$clust.distr, round.pos))
  cat("\n")
  
    
  if(!is.null(object$e.values))
  {
    cat("Eigen values:\n")
    print(round(object$e.values$values, round.pos))
    cat("\n")
  }
  
  cat("Measures with max.loading in factors:\n")
  print(names(object$factors)[-1])
  cat("\n")
  
  cat("Summary of measures with max.loading in factors:\n")
  summary.fact = (round(apply(object$factor[,-1], 2, summary), round.pos))
  std = round(apply(object$factor[,-1], 2, sd), round.pos)
  print(rbind(summary.fact, std))
  cat("\n")
  
  cat("Summary of measures with max.loading in factors by cluster:\n")
  cat("\n")

  clust.num = length(object$clust.distr)
  for(i_clust in 1:clust.num)
  {
    clust.id = object$clusters$ID[object$clusters$cluster == i_clust]
    clust.fact = object$factors[which(object$factors$output %in% clust.id),]
    cat("Cluster", i_clust, "\n")
    summary.clust.fact = round(apply(clust.fact[,-1], 2, summary), round.pos)
    std = round(apply(clust.fact[-1], 2, sd),round.pos)
    print(rbind(summary.clust.fact,std))
    cat("\n")
  }
  cat('\nIf you report these results, please cite:\nSylvestre MP, et al. (2006). Classification of patterns of delirium severity scores over time in an elderly population. 
International Psychogeriatrics,18(4), 667-680. doi:10.1017/S1041610206003334.\n')

}
