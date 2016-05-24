compactness <- function(g, dist=NULL, wien=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist=distanceMatrix(g)
  }
  if(is.null(wien)){
    wien <- wiener(g,dist)
  }
  nV <- numNodes(g)
  return (4*wien/(nV*(nV-1)))
}
