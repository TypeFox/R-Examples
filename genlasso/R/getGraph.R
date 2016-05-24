getGraph <- function(D) {
  mat = abs(crossprod(D))
  diag(mat) = 0
  return(graph.adjacency(mat, mode=c("upper")))
}
