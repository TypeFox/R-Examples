dendro.variables <- function(data, dist.variables.method=c("associationMeasures", "ClustOfVar"), associationFun=association, check.psd=TRUE){
  
  dist.variables.method <- match.arg(dist.variables.method)
  if(dist.variables.method == "associationMeasures"){
    S <- similarity.variables(data, associationFun=associationFun, check.psd=check.psd)
    D.variables <- as.dist(sqrt(1 - S))
    dend <- as.dendrogram(hclust(D.variables))
  }

  else if(dist.variables.method == "ClustOfVar"){
     dc <- sapply(data, data.class)
    if(any(dc == "numeric"))
      X.quanti <- data[,dc == "numeric"]
    else
      X.quanti <- NULL
    if(all(dc == "numeric"))
      X.quali <- NULL
    else
      X.quali <- data[,dc != "numeric"]
    dend <- as.dendrogram(ClustOfVar::hclustvar(X.quanti, X.quali))
  }
  
  return(dend)
}
