# This function constructs the graph Laplacian matrices from adjacency matrices
# Author: Sen Zhao
# Email: sendavid7@gmail.com

make.L <- function(adj, normalize.Laplacian = FALSE){
  if(isSymmetric(adj) == FALSE){
    stop("Error: The adjacency matrix needs to be symmetric.")
  }
  if(max(abs(adj)) > 1){
    stop("Error: The weight in the adjacency matrix needs to be between -1 and 1.")
  }
  L <- -adj
  diag(L) <- 0
  diag(L) <- -rowSums(L)
  if(normalize.Laplacian){
    diag(L)[diag(L) == 0] <- 1
    L <- diag(1 / sqrt(diag(L))) %*% L %*% diag(1 / sqrt(diag(L)))
  }
  return(L)
}
