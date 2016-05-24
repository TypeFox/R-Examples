`my.positive.definite.solve` <-
function(A, eps = 1e-15 ) {
  h <- eigen(A,symmetric=TRUE)
  ind <-  (1:(dim(A))[1])[h$values >= eps]
  return( ( h$vectors[,ind] %*% diag(1/h$values[ind],length(ind)) %*% t( h$vectors[,ind])))
}

