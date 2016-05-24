geometricArithmetic3 <- function(g, dist=NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  M <- adjacencyMatrix(g)

  sum(.edgeApply(M, function(i, j) {
    to_i <- .edgeApply(M, function(k, l) {
      min(dist[k, i], dist[l, i])
    }, dupls=FALSE)
    to_j <- .edgeApply(M, function(k, l) {
      min(dist[k, j], dist[l, j])
    }, dupls=FALSE)
    m_i <- sum(to_i < to_j)
    m_j <- sum(to_i > to_j)
    2 * sqrt(m_i * m_j) / (m_i + m_j)
  }, dupls=FALSE))
}
