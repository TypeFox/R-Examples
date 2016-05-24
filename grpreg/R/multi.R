multiX <- function(X, m) {
  p <- ncol(X)
  n <- nrow(X)
  A <- matrix(0, m*n, m*p)
  for (i in 1:m) {
    A[m*(1:n)-i+1, m*(1:p)-i+1] <- X
  }
  cbind(matrix(as.numeric(diag(m)),m*n,m,byrow=TRUE)[,2:m],A)
}
multiY <- function(y) {
  as.numeric(t(y))
}
