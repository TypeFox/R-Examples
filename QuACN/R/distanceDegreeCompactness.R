distanceDegreeCompactness <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  if (is.null(dist))
    dist <- distanceMatrix(g)

  W <- wiener(g, dist=dist)
  distdeg <- rowSums(dist)
  ecc <- apply(dist, 1, max)
  center <- which.min(ecc)
  radius <- ecc[[center]]

  qks <- sapply(1:radius, function(k) {
    vs <- which(dist[center,] == k)
    sum(distdeg[vs])
  })

  2*W * log2(2*W) - sum(qks * log2(qks))
}
