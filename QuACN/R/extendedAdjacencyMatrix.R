extendedAdjacencyMatrix <- function(g) {
  stopifnot(.validateGraph(g))
  deg.vec <- graph::degree(g)
  n_nodes <- length(deg.vec)
  ExtAdjMat <- matrix(0, nrow=n_nodes, ncol=n_nodes, byrow=TRUE)
  adj.mat <- adjacencyMatrix(g) 
  for(i in 1:n_nodes) {
    for(j in 1:n_nodes) {
      if(adj.mat[i,j] != 0) {
        ExtAdjMat[i,j] <- ((deg.vec[i]/deg.vec[j]) + (deg.vec[j]/deg.vec[i]))/2
      }
    }
  }
  ExtAdjMat
}
