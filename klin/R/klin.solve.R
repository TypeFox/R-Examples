`klin.solve` <-
function(A,b) {
  ## dimensions of matrices
  if (!is.list(A))
    stop("first argument should be a list")
  if (!all(sapply(A,function(m) { is(m,"Matrix") || is(m,"matrix") })))
    stop("first argument should be a list of matrices (matrix or Matrix)")
  m <- sapply(A,nrow)
  n <- sapply(A,ncol)
  K <- length(A)
  ## check consistency
  if (!is.vector(b))
    stop("second argument should be a vector")
  if (prod(m)!=length(b))
    stop("incompatible dimensions")
  if (any(n!=m))
    stop("matrices need to be square")
  ## initialize
  B <- b
  ## loop
  for (k in K:1) {
    B <- Matrix(as.vector(B),m[k],prod(c(m[incseq(1,k-1)],n[incseq(k+1,K)])))
    B <- solve(A[[k]],B)
    B <- t(B)
  }
  as.vector(B)
}
