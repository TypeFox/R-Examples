distanceDegreeEquality <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  if (is.null(dist))
    dist <- distanceMatrix(g)

  distdeg <- rowSums(dist)
  part <- table(distdeg)
  pis <- part / numNodes(g)

  -sum(pis * log2(pis))
}
