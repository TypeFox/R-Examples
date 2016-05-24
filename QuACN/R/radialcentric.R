radialCentric <- function(g, dist=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  ecc <- apply(dist,1,max)
  pis <- table(ecc)/numNodes(g)

  return ((-1)*sum(pis*log2(pis)))
}
