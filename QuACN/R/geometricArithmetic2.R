geometricArithmetic2 <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  sum(.edgeApply(g, function(i, j) {
    to_i <- dist[, i]
    to_j <- dist[, j]
    n_i <- sum(to_i < to_j)
    n_j <- sum(to_i > to_j)
    2 * sqrt(n_i * n_j) / (n_i + n_j)
  }, dupls=FALSE))
}
