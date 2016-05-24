konstantinova <- function(g, dist=NULL){
   # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist = distanceMatrix(g)
  }
  V = numNodes(g)
  d <- rowSums(dist)
  p <- unlist(dist/d)
  p <- p[p!=0]
  
  return(-sum(p*log2(p)))
}
