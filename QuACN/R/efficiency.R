efficiency <- function(g, dist=NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  n <- numNodes(g)
  m_Clique <- choose(n, 2)
  sum_path <- sum(sapply(1:(n-1), function(j) (n-j)/j))
  E_path <- 1/m_Clique * sum_path

  inv_dist <- 1/dist
  inv_dist[lower.tri(inv_dist, TRUE)] <- 0
  sumCE <- sum(inv_dist)

  E <- 1/m_Clique * sumCE
  (4*(E - E_path)*(1-E)) / ((1-E_path)^2)
}
