bonchev3 <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))

  if (is.null(dist))
    dist <- distanceMatrix(g)

  rho <- max(dist)
  ki <- table(dist)[2:(rho+1)]
  nV <- numNodes(g)
  pis <- ki / nV / (nV - 1)

  -sum(pis * log2(pis))
}
