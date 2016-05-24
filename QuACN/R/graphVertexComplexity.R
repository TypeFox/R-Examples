graphVertexComplexity <- function(g, dist=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist=distanceMatrix(g)
  }
  #names of nodes
  nam <- nodes(g)
  #number of vertices
  nV <- numNodes(g)
  #determine distance distributen for each vertex
  distd <- lapply(nam,function(n){
    table(dist[n,])#,exclude=0)
  })
  names(distd) <- nam
  vic <- sapply(distd,function(dd){
    pis <- dd/nV
    -sum(pis*log2(pis))
  })
  return(sum(vic)/nV)
}
