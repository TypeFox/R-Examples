estrada <- function(g) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  n <- numNodes(g)
  m <- numEdges(g)

  M <- adjacencyMatrix(g)
  EV <- as.double(eigen(M, only.values=TRUE)$values)
  sum(exp(EV))
}
