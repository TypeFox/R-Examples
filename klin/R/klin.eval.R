`klin.eval` <-
function(A, x, transpose=FALSE) {
  ## dimensions of matrices
  m <- sapply(A,nrow)
  n <- sapply(A,ncol)
  K <- length(A)
  ## degenerate case
  if (K==0)
    return(x)
  ## switch n and m if we are working with t(A)
  if (transpose) {
    tmp <- m
    m <- n
    n <- tmp
  }
  ## check consistency
  if (!is.list(A))
    stop("first argument should be a list")
  if (!is.vector(x))
    stop("second argument should be a vector")
  if (prod(n)!=length(x))
    stop("incompatible dimensions")
  ## initialize
  X <- x
  ## loop
  for (k in K:1) {
    X <- Matrix(as.vector(X),n[k],prod(c(n[incseq(1,k-1)],m[incseq(k+1,K)])))
    if (transpose)
      Y <- crossprod(A[[k]],X)
    else
      Y <- A[[k]] %*% X
    X <- t(Y)
  }
  as.vector(X)
}

