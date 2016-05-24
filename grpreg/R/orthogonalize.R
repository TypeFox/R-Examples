orthogonalize <- function(X, group) {
  n <- nrow(X)
  J <- max(group)
  T <- vector("list", J)
  XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
  XX[,which(group==0)] <- X[,which(group==0)]
  for (j in seq_along(numeric(J))) {
    ind <- which(group==j)
    if (length(ind)==0) next
    SVD <- svd(X[, ind, drop=FALSE], nu=0)
    r <- which(SVD$d > 1e-10)
    T[[j]] <- sweep(SVD$v[,r,drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
    XX[,ind[r]] <- X[,ind]%*%T[[j]]
  }
  nz <- !apply(XX==0,2,all)
  XX <- XX[, nz, drop=FALSE]
  attr(XX, "T") <- T
  attr(XX, "group") <- group[nz]
  XX
}
unorthogonalize <- function(b, XX, group) {
  ind <- !sapply(attr(XX, "T"), is.null)
  T <- bdiag(attr(XX, "T")[ind])
  ind0 <- c(1, 1+which(group==0))
  rbind(b[ind0,,drop=FALSE], as.matrix(T %*% b[-ind0,,drop=FALSE]))
}
