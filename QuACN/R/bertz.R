bertz <- function(g,dist=NULL,deg=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))

  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  if(is.null(deg)){
    deg <- graph::degree(g)
  }
  Ni <- .cardNi(g,dist,deg)
  N <- sum(Ni)
  return((2*N*log2(N))- sum(Ni*log2(Ni)))
}
