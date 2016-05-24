GenerateCliquesCovariance <- function(ncliques, cliquesize, theta) {
  # Generates a block diagonal positive definite matrix with randomly signed
  # non-zero elements and condition number equal to p=ncliques*cliquesize
  #
  # Args:
  #   ncliques: number of cliques
  #   cliquesize: size of clique
  #   theta: magnitude of non-zeros
  p <- ncliques * cliquesize
  sizes <- rep(cliquesize, ncliques)
  Sigma <- matrix(0, p, p)
  lims <- c(0, cumsum(sizes))
  for (i in seq(ncliques)) {
    ii <- (lims[i] + 1):lims[i+1]
    signs <- 2 * (rnorm(sizes[i] * (sizes[i] - 1) / 2) < 0) - 1
    Sigma[ii, ii][upper.tri(Sigma[ii, ii])] <- signs * theta
  }
  Sigma <- Sigma + t(Sigma)
  eig <- eigen(Sigma, symmetric=T)
  shift <- (max(eig$val) - p * min(eig$val)) / (p - 1)
  cat("Shifting eigenvalues by ", shift, fill=T)
  diag(Sigma) <- diag(Sigma) + shift

  # compute symmetric square root for generating data
  A <- eig$vect %*% diag(sqrt(eig$val + shift)) %*% t(eig$vect)
  list(Sigma=Sigma, A=A, shift=shift)
}
