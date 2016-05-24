balabanJ <- function(g, dist=NULL){
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))

  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  
  sum.dist <- rowSums(dist)
  edges <- .edgesNoDupls(g)

  Ji<-sapply(1:length(edges), function(i){
    start <- names(edges[i])
    targets <- names(edges[[i]])
    1/sqrt(sum.dist[start] * sum.dist[targets])
  })
  
  nE <- numEdges(g)
  nV <- numNodes(g)

  return ( (nE/(nE-nV+2)) * sum(unlist(Ji)))
}
