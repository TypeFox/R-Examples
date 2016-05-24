cluster<- function(x,clusters)
{
  if ("Simpsons" %in% class(x)) stop("Simpsons object needed")
  
  Nclusters <- x$Nclusters
  
  splList <- split(x$alldata,x$alldata$clusterid)
  names(splList) <- paste("Cluster",1:Nclusters)
  
  if (!missing(clusters))
  {
    splList <- splList[clusters]
    mat <- do.call(rbind,splList)
    rownames(mat) <- NULL
    return(mat)
  }
  
  return(splList)
}

