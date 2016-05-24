graphDistanceComplexity <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  distdeg <- rowSums(dist)

  lis <- sapply(nodes(g), function(i) {
    ecc <- max(dist[i,])
    -sum(sapply(1:ecc, function(j) {
      aij <- sum(dist[i,] == j)
      p <- j / distdeg[[i]]
      aij * p * log2(p)
    }))
  })

  sum(lis) / length(lis)
}
