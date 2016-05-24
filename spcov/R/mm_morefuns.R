# "Internal" functions for generalized gradient descent

g <- function(Sigma, S, Sigma0, inv.Sigma0=NULL, log.det.Sigma0=NULL) {
  # the differentiable part of the objective
  p <- nrow(Sigma)
  ind.diag <- 1 + 0:(p - 1) * (p + 1)
  if (is.null(log.det.Sigma0))
      log.det.Sigma0 <- LogDet(Sigma0)
  if (is.null(inv.Sigma0))
    log.det.Sigma0 + sum((solve(Sigma0, Sigma) + solve(Sigma, S))[ind.diag]) - p
  else
    log.det.Sigma0 + sum((inv.Sigma0 %*% Sigma + solve(Sigma, S))[ind.diag]) - p
}

ComputeGradientOfg <- function(Sigma, S, Sigma0, inv.Sigma0=NULL) {
  # computes grad g
  inv.Sigma <- solve(Sigma)
  if (is.null(inv.Sigma0))
    solve(Sigma0) - inv.Sigma %*% S %*% inv.Sigma
  else
    inv.Sigma0 - inv.Sigma %*% S %*% inv.Sigma
  # note: this can be sped up by using solve(Sigma, X), since S=t(X)%*%X/n
  # and if S=t(X)%*%X/n + eps*diag(p), we can use Sherman-Morrison
}

ComputeLikelihood <- function(Sigma, S) {
    # (1/n times) the log-likelihood of a multivariate normal
  p <- nrow(Sigma)
  ind.diag <- 1 + 0:(p - 1) * (p + 1)
  -(1 / 2) * (LogDet(Sigma) + sum(solve(Sigma, S)[ind.diag]))
}

ComputePenalty <- function(Sigma, lambda) {
  # the non-differentiable part of the objective
  # Args:
  #  lambda: penalty parameter.  Can be a scalar or a matrix of the same
  # dimensions as Sigma
  if (length(lambda) != 1 & any(dim(lambda) != dim(Sigma)))
    stop("lambda must be a scalar or else a square matrix of size nrow(Sigma).")
  sum(lambda * abs(Sigma))
}

ComputeMajorizerObjective <- function(Sigma, S, Sigma0, lambda,
                                      inv.Sigma0=NULL, log.det.Sigma0=NULL) {
  # the majorizer's objective function
  g(Sigma, S, Sigma0, inv.Sigma0=inv.Sigma0, log.det.Sigma0=log.det.Sigma0) + ComputePenalty(Sigma, lambda)
}

LogDet <- function(M) {
  # note: determinant() is better than log(det()) when det is very close to zero
  log.det <- determinant(M, logarithm=TRUE)
  if (log.det$sign == -1)
    return(NA)
  as.numeric(log.det$mod)
  #log(det(M))
}


