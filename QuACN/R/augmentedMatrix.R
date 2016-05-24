augmentedMatrix <- function(g) {
  stopifnot(.validateGraph(g))

  dist.mat <- distanceMatrix(g)
  n_nodes <- nrow(dist.mat)
  Aug_Mat <- matrix(0, nrow=n_nodes, ncol=n_nodes, byrow=TRUE)
  deg.vec <- graph::degree(g)
  for(i in 1:n_nodes) {
    for(j in 1:n_nodes) {
      if(i != j) {
        Aug_Mat[i,j] <- deg.vec[j]/(2^{dist.mat[i,j]})
      }
      else {
        Aug_Mat[i,i] <- deg.vec[i]
      }
    }
  }
  Aug_Mat
}
