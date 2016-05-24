distanceDegreeCentric <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  n <- numNodes(g)

  sum.dist <- rowSums(dist)
  max.dist <- apply(dist, 1, max)

  freq <- table(sum.dist, max.dist)
  freq <- freq[freq != 0]
  p <- freq / n

  -sum(p * log2(p))
}
