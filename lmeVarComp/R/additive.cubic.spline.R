additive.cubic.spline <-
function(x, nbasis, weights = NULL)
{
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }
  
  p <- ncol(x)
  interval <- c(0, 1)
  knots <- seq(0, 1, length.out = nbasis - 2L)[2L : (nbasis - 3L)]
  
  Z <- lapply(seq_len(p), function(j) 
    cubic.spline(x[, j], knots, interval))
  S <- lapply(seq_len(p), function(j) 
    cubic.spline.penalty(knots, interval, 2L))
  U <- lapply(seq_len(p), function(j)
    mnls(Z[[j]], cbind(1, x[, j])))
  
  Z <- do.call(cbind, Z)
  S <- do.call(block.diag, S)
  U <- do.call(block.diag, U)
  U[is.na(U)] <- 0
  V <- crossprod(Z, Z * (weights ^ 2)) %*% U
  V <- qr.Q(qr(V, LAPACK = TRUE), TRUE)[, -seq_len(ncol(U))]
    # then ZU and ZV are column orthogonal
  
  X <- cbind(1, x) * weights  # spans the same space as ZU does
  Z <- Z %*% V * weights      # satsifies X^T Z = 0
  S <- crossprod(V, S %*% V)
  S <- (S + t(S)) / 2

  list(X = X, Z = Z, S = S)   # X == cbind(1, x) * weights
}
