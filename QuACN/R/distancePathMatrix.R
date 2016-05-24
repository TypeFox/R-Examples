distancePathMatrix <- function(g) {
  stopifnot(.validateGraph(g))
  
  dist.mat <- distanceMatrix(g)
  n_nodes <- nrow(dist.mat)
  Dis_Path_Mat <- matrix(0, nrow=n_nodes, ncol=n_nodes, byrow=TRUE)
  for(i in 1:n_nodes) {
    for(j in 1:n_nodes) {
      Dis_Path_Mat[i,j] <- choose(dist.mat[i,j] + 1, 2)
    }
  }
  Dis_Path_Mat
}
