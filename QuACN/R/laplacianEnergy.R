laplacianEnergy <- function(g) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  n <- numNodes(g)
  m <- numEdges(g)

  lap <- laplaceMatrix(g)
  EV <- as.double(eigen(lap, only.values=TRUE)$values)
  sum(abs(EV - 2*m/n))
}
