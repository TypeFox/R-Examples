bonchev1 <- function(g,dist=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  rho <- max(dist)
  ki <- table(dist)[2:(rho+1)]
  nV <- numNodes(g)
  pis <- 2*ki/nV^2
  In <- (pis*log2(pis))
  return (((-1)/nV)*log2(1/nV)-sum(In))
}
