laplacianEstrada <- function(g) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")

  n <- numNodes(g)
  m <- numEdges(g)

  lap <- laplaceMatrix(g)
  EV <- as.double(eigen(lap, only.values=TRUE)$values)
  sum(exp(EV))
}
