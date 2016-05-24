dobrynin<- function(g, dist=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }

  nN <- numNodes(g)
  nam <- nodes(g)
  res <- list()
  res$eccentricityVertex <- apply(dist,1,max)
  names(res$eccentricityVertex) <- nam
  res$eccentricityGraph <- sum(res$eccentricityVertex)
  res$avgeccOfG <- res$eccentricityGraph/nN
  res$ecentricVertex <- abs(res$eccentricityVertex - res$avgeccOfG)
  names(res$ecentricVertex) <- nam
  res$ecentricGraph <- sum(res$ecentricVertex/nN)
  res$vertexCentrality <- rowSums(dist)
  names(res$vertexCentrality) <- nam
  res$graphIntegration <- sum(res$vertexCentrality)/2
  res$unipolarity <- min(res$vertexCentrality)
  res$vertexDeviation <- res$vertexCentrality - res$unipolarity
  names(res$vertexDeviation) <- nam
  res$variation <- max(res$vertexDeviation)
  res$centralization <- sum(res$vertexDeviation)
  res$avgDistance <- 2*res$graphIntegration/nN
  res$distVertexDeviation <- abs(res$vertexCentrality - res$avgDistance)
  names(res$distVertexDeviation) <- nam
  res$meanDistVertexDeviation <- mean(res$distVertexDeviation)
  
  return(res)
}  
