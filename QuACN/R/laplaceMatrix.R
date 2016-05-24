laplaceMatrix <- function(g) {
  stopifnot(.validateGraph(g))
  adj.mat <- adjacencyMatrix(g)
  D <- diag(rowSums(adj.mat, na.rm = FALSE, dims = 1))
  Lap_Mat <- D - adj.mat
  Lap_Mat
}
