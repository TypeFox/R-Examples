find_log_integrand<- function(pars, modfr)
{
  n_fixed <- ncol(modfr$X)
  n_random <- length(modfr$reTrms$theta)

  Zt <- modfr$reTrms$Zt
  LambdatTheta <- modfr$reTrms$Lambdat
  theta <- pars[1:n_random]

  LambdatTheta@x[] <- theta[modfr$reTrms$Lind]
  LambdatThetaZt <- LambdatTheta %*% Zt

  family <- modfr$family
  y <- model.response(modfr$fr)
  weights <- model.weights(modfr$fr)
  if(length(weights) == 0)
    weights <- rep(1, length(y))

  beta <- pars[(n_random + 1) : (n_random + n_fixed)]
  XBeta <- as.numeric(modfr$X %*% beta)

  log_integrand <- function(u) {
    utLambdatThetaZt <- Matrix::crossprod(u, LambdatThetaZt)

    utLambdatThetaZt_num <- as.numeric(utLambdatThetaZt)
    eta_num <- rep(XBeta, each = ncol(u)) + utLambdatThetaZt_num
    mu_num <- family$linkinv(eta_num)
    y_rep <- rep(y, each = ncol(u))
    weights_rep <- rep(weights, each = ncol(u))
    dev_resids_num <- family$dev.resids(y_rep, mu_num, weights_rep)
    dev_resids <- matrix(dev_resids_num, nrow = ncol(u))
    - rowSums(dev_resids) / 2 + colSums(dnorm(u, log = TRUE))
  }
  log_integrand
}

transform_normal <- function(z, normal_approx) {
  C <- solve(chol(normal_approx$precision))
  normal_approx$mean + crossprod(C, z)
}

find_log_weights <- function(z, log_integrand, normal_approx) {
  u <- transform_normal(z, normal_approx)
  log_det_prec <- determinant(as.matrix(normal_approx$precision), logarithm = TRUE)$modulus
  normal_part <- colSums(dnorm(z, log = TRUE)) + log_det_prec/2
  log_weights <- log_integrand(u) - normal_part
  log_weights
}

find_lfun_IS <- function(modfr, devfun_laplace, nIS) {
  n <- nrow(modfr$reTrms$Zt)
  z_poss <- matrix(rnorm(n * nIS), nrow = n)

  lfun <- function(pars) {
    normal_approx <- compute_normal_approx(pars, devfun_laplace)
    log_integrand <- find_log_integrand(pars, modfr)
    log_weights <- find_log_weights(z_poss, log_integrand, normal_approx)
    b <- mean(log_weights)
    log(mean(exp(log_weights - b))) + b
  }

  lfun
}
