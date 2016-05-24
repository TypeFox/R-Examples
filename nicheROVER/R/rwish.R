rwish <-
function(n, Psi, nu, inv = FALSE) {
  if(inv) Psi <- solve(Psi)
  U <- chol(Psi)
  d <- nrow(Psi)
  ans <- array(0, dim = c(d, d, n))
  if(!is.null(dimnames(Psi))) dimnames(ans) <- c(dimnames(Psi), list(NULL))
  ans[rep(upper.tri(Psi), n)] <- rnorm(n*d*(d-1)/2)
  ans[rep(!lower.tri(Psi, diag = FALSE) &
            !upper.tri(Psi, diag = FALSE), n)] <- sqrt(rchisq(n*d, df = nu-1:d+1))
  for(ii in 1:n) {
    tmp <- ans[,,ii] %*% U
    if(inv) tmp <- backsolve(tmp, diag(d), transpose = TRUE)
    ans[,,ii] <- crossprod(tmp)
  }
  ans
}
