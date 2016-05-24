niiw.post <-
function(nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, burn) {
  # sufficient statistics
  d <- ncol(X)
  N <- nrow(X)
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # local variables
  mu <- rep(0, d)
  Sigma <- matrix(0, d, d)
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- floor(nsamples*burn)
  mu <- mu0
  # output
  mu.out <- matrix(NA, nsamples, d)
  Sigma.out <- array(NA, dim = c(d, d, nsamples))
  # main for loop
  for(ii in (-burn+1):nsamples) {
    # sample Sigma
    Psi2 <- N * ((mu-Xbar) %o% (mu-Xbar)) + S + Psi
    nu2 <- N+nu
    Sigma <- matrix(rwish(1, Psi2, nu2, inv = TRUE), d, d)
    # sample mu
    Sigma2 <- Sigma/N
    if(!all(Omega == 0)) {
      B <- Sigma2 %*% solve(Omega + Sigma2)
    } else B <- matrix(0, d, d)
    IB <- diag(d)-B
    lambda2 <- c(IB %*% Xbar)
    if(!all(Omega == 0)) lambda2 <- lambda2 + c(B %*% lambda)
    Omega2 <- IB %*% Sigma2
    mu <- c(rmvnorm(1, lambda2, Omega2))
    # store
    if(ii > 0) {
      mu.out[ii,] <- mu
      Sigma.out[,,ii] <- Sigma
    }
  }
  list(mu = mu.out, Sigma = Sigma.out)
}
