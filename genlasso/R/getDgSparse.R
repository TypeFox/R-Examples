getDgSparse <- function(graph) {
  E = get.edgelist(graph)
  A = rbind(cbind(Seq(1,nrow(E)),E[,1],-1),cbind(Seq(1,nrow(E)),E[,2],1))
  D = sparseMatrix(i=A[,1], j=A[,2], x=A[,3], dims=c(ecount(graph),vcount(graph)))
  return(D)
}
