###
# A collection of helper functions for testing, help create random examples.
###

rConnectedAdjMatrix = function(n,p) {
  weights = runif(n*(n-1)/2)
  g = minimum.spanning.tree(graph.full(n), weights=weights)
  adjMatrix = as.matrix(get.adjacency(g))
  adjMatrix = (upper.tri(matrix(0,n,n)) &
                 matrix(sample(c(T, F), n^2,
                               replace=T, prob=c(p, 1-p)), ncol=n)) | adjMatrix
  adjMatrix = 1*(adjMatrix | t(adjMatrix))
  return(adjMatrix)
}

rDirectedAdjMatrix = function(n,p) {
  return(1*(upper.tri(matrix(0,n,n)) & matrix(sample(c(T, F), n^2, replace=T,
                                                     prob=c(p, 1-p)), ncol=n)))
}