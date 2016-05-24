VertexSimilarity <- function(m,alpha=0.97){
  requireNamespace("igraph")
  network=igraph::graph.adjacency(m,mode="undirected")
  n <- igraph::vcount(network)
  mat <- matrix(0,nrow = n,ncol = n)
  id <- matrix(0,nrow=n,ncol=n)
  sim <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    mat[i,i] <- igraph::degree(network,i,mode="all")
    id[i,i] <- 1
    sim[i,i] <- 1
  }
  mat <- solve(mat)
  eig <- max(unlist(eigen(m,only.values = TRUE)))
  e <- igraph::ecount(network)
  sim <- mat%*%solve((id-(alpha/eig)*m))%*%mat
  sim <- 2*e*eig*sim
  return(sim)
}
