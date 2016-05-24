convertGraph <- function(Graph) {
    # FUNCTION
  #   Convert graphs into adjacency matrices
  # INPUT
  #   Graph:  array p x max(dmax) x length(K)
  # (dmax, K :see selectFast)
  # When length(K)=1, a matrix
  # Graph[a,,iK]: indexes of the nodes connected to the node a,
  #   for the value iK of K;
  # Graph[a,1,iK]=0 if there is no node connected to the node a.
  # OUTPUT
  #   AdjG: array p x p x length(K)
  #  AdjG[,,iK] is a symmetric matrix, with diagonal equal to zero
   # When length(K)=1, a matrix
  # Note:
  # Chapeau a convertGraph, pour tenir compte du cas length(K)=1
  # quand appele par l'utilisateur
  # CALLED BY
  #   End-user function.
  # ---------------------------------------------------------------
  # Si length(K)=1, Graph est une matrice
  if (is.matrix(Graph)) {
    Graph <- array(Graph, c(dim(Graph)[1], dim(Graph)[2],1),
                   dimnames=list(dimnames(Graph)[[1]],NULL,NULL))
  }
  res <- convGraph(Graph)
if (dim(res)[[3]]  ==1) {
  return(as.matrix(res[,,1]))
}
  else
    return(res)
} # fin convertGraph
  
convGraph <- function(Graph) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Convert graphs into adjacency matrices
  # INPUT
  #   Graph:  array p x max(dmax) x length(K)
  # (dmax, K :see selectFast)
  # Graph[a,,iK]: indexes of the nodes connected to the node a,
  #   for the value iK of K;
  # Graph[a,1,iK]=0 if there is no node connected to the node a.
  # OUTPUT
  #   AdjG: array p x p x length(K)
  #  AdjG[,,iK] is a symmetric matrix, with diagonal equal to zero
  # CALLED BY
  #   calcLarsNEW
  # ---------------------------------------------------------------
  lK <- dim(Graph)[3]
  p <- dim(Graph)[1]
  AdjG <- array(0,c(p, p, lK))
  for (iK in 1:lK) {
    for (a in 1:p) {
      ind <- Graph[a,, iK]!=0
      if (sum(ind)>0) {
        AdjG[a, Graph[a, ind,  iK], iK] <- 1
      }
    }
  }
  dimnames(AdjG) <- list(dimnames(Graph)[[1]],
                         dimnames(Graph)[[1]],
                         dimnames(Graph)[[3]])
  return(AdjG)
}
