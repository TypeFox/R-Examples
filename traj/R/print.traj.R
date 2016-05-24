print.traj <-
function(x, round.pos = 2, ...)
{
  cat("Number of observations: ", nrow(x$data),"\n\n")
  
  cat("Cluster distribution:\n")
  print(round(x$clust.distr, round.pos))
  cat("\n")
  
  cat("Measures with max.loading in factors: ", names(x$factors)[-1])
  cat("\n")
  cat('\nIf you report these results, please cite:\nSylvestre MP, et al. (2006). Classification of patterns of delirium severity scores over time in an elderly population. 
International Psychogeriatrics,18(4), 667-680. doi:10.1017/S1041610206003334.\n')
}
