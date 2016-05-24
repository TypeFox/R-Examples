eigenl <-
function(pm)
{
  #"pm" is a matrix
  #If matrix is quadratic (i.e. number of age classes > 1), calculate
  #eigenvalues and vectors the normal way
  if(dim(pm)[1] == dim(pm)[2])
  {
    #The dominant eigenvalue
    ret <- list(lambda = as.numeric(eigen(pm, only.values = TRUE)$values[1]))
    #The right eigenvector of pm
    ret$u <- abs(eigen(pm)$vectors[, 1])
    #The left eigenvector of pm is found by taking the right egenvector of
    #the transposed pm
    ret$v <- abs(eigen(t(pm))$vectors[, 1])
    #Scaling
    ret$u <- ret$u / sum(ret$u)
    ret$v <- ret$v / sum(ret$u * ret$v)
  }
  #Otherwise, define lambda, and the stable age distribution (u) and
  #reproductive values (v) as below
  else
  {
    ret <- list(lambda = colSums(pm))
    ret$u <- 1
    ret$v <- 1
  }
  #Output
  ret
}
