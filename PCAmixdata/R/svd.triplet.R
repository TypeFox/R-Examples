svd.triplet<-function (X, row.w = NULL, col.w = NULL, ncp = Inf) 
{
  if (is.null(row.w)) 
    row.w <- rep(1/nrow(X), nrow(X))
  if (is.null(col.w)) 
    col.w <- rep(1, ncol(X))
  ncp <- min(ncp, nrow(X) - 1, ncol(X))
  row.w = row.w/sum(row.w)
  X = sweep(X, 2, sqrt(col.w), FUN = "*")
  X = sweep(X, 1, sqrt(row.w), FUN = "*")
  if (ncol(X) < nrow(X)) {
    svd.usuelle <- svd(X)
    U <- svd.usuelle$u[, 1:ncp, drop = FALSE]
    V <- svd.usuelle$v[, 1:ncp, drop = FALSE]
    if (ncp > 1) {
      mult <- sign(apply(V, 2, sum))
      mult[mult == 0] <- 1
      U <- sweep(U, 2, mult, FUN = "*")
      V <- sweep(V, 2, mult, FUN = "*")
    }
    U <- sweep(as.matrix(U), 1, sqrt(row.w), FUN = "/")
    V <- sweep(as.matrix(V), 1, sqrt(col.w), FUN = "/")
  }
  else {
    svd.usuelle <- svd(t(X))
    U <- svd.usuelle$v[, 1:ncp, drop = FALSE]
    V <- svd.usuelle$u[, 1:ncp, drop = FALSE]
    mult <- sign(apply(V, 2, sum))
    mult[mult == 0] <- 1
    V <- sweep(V, 2, mult, FUN = "*")
    U <- sweep(U, 2, mult, FUN = "*")
    U <- sweep(U, 1, sqrt(row.w), FUN = "/")
    V <- sweep(V, 1, sqrt(col.w), FUN = "/")
  }
  vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
  num <- which(vs[1:ncp] < 1e-15)
  if (length(num) == 1) {
    U[, num] <- U[, num] * vs[num]
    V[, num] <- V[, num] * vs[num]
  }
  if (length(num) > 1) {
    U[, num] <- sweep(U[, num], 2, vs[num], FUN = "*")
    V[, num] <- sweep(V[, num], 2, vs[num], FUN = "*")
  }
  res <- list(vs = vs, U = U, V = V)
  return(res)
}