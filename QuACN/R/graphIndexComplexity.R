graphIndexComplexity <- function(g) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  n <- numNodes(g)
  M <- adjacencyMatrix(g)
  EV <- as.double(eigen(M, only.values=TRUE)$values)
  r_max <- max(EV)
  cr <- (r_max - 2*cos(pi/(n+1))) / ((n-1) - 2*cos(pi/(n+1)))

  4 * cr * (1-cr)
}
