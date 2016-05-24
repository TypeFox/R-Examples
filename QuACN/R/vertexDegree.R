vertexDegree <- function(g, deg=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(deg)){
    deg <- graph::degree(g)
  }
  Nikv <- table(deg)
  pis <- Nikv/numNodes(g)
  return ((-1)*sum(pis*log2(pis)))
}
