clusterAnalysis <-
function(stmobj, labels_number=3){
  labels <- stm::labelTopics(stmobj)
  theta <- stmobj$theta
  d <- stats::dist(stats::cor(theta))
  clust <- stats::hclust(d, method="complete")
  K <- stmobj$settings$dim$K
  clust$labels <- rep(NA,K)
  for(i in seq(K)){
    l <- labels$frex[i,]
    l <- paste(l[1:labels_number], collapse=', ')
    clust$labels[i] <- l
  }
  return(clust)
}
