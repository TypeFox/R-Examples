localClusteringCoeff <- function(g, deg=NULL){
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(deg)){ 
    deg <- graph::degree(g)
  }
  #lcc <- (2*length(edges(g)))/(deg*(deg-1))
  lcc <- igraph::transitivity(.G2IG(g), type='local')
  names(lcc) <- nodes(g)
  lcc[is.na(lcc)] <- 0

  return(lcc)
}
