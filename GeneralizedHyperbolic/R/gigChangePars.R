### Change parameterizations of the generalized inverse Gaussian distribution
gigChangePars <- function (from, to, param, noNames = FALSE) {

  param <- as.numeric(param)

  # If lambda is ommitted from the param vector, it defaults to 1
  if (length(param) == 2)
    param <- c(param, 1)

  if (length(param) != 3)
    stop("param vector must contain 3 values")

  if (! from %in% 1:4)
    stop("the argument 'from' must be either 1, 2, 3 or 4")

  if (! to %in% 1:4)
    stop("the argument 'to' must be either 1, 2, 3 or 4")

  param <- as.numeric(param)

  lambda <- param[3]

  if (from == 1) {
    chi <- param[1]
    psi <- param[2]

    if (chi <= 0)
      stop("chi must be greater than zero")

    if (psi <= 0)
      stop("psi must be greater than zero")
  }

  if (from == 2) {
    delta <- param[1]
    gamma <- param[2]

    if (delta <= 0)
      stop("delta must be greater than zero")

    if (gamma <= 0)
      stop("gamma must be greater than zero")
  }

  if (from == 3) {
    alpha <- param[1]
    beta <- param[2]

    if (alpha <= 0)
      stop("alpha must be greater than zero")

    if (beta <= 0)
      stop("beta must be greater than zero")
  }

  if (from == 4) {
    omega <- param[1]
    eta <- param[2]

    if (omega <= 0)
      stop("omega must be greater than zero")

    if (eta <= 0)
      stop("eta must be greater than zero")
  }

  if (from == 1 && to == 2) {
    delta <- sqrt(chi)
    gamma <- sqrt(psi)
    output <- c(delta = delta, gamma = gamma, lambda = lambda)
  }

  if (from == 1 && to == 3) {
    alpha <- sqrt(psi / chi)
    beta <- sqrt(chi * psi)
    output <- c(alpha = alpha, beta = beta, lambda = lambda)
  }

  if (from == 1 && to == 4) {
    omega <- sqrt(chi * psi)
    eta <- sqrt(chi / psi)
    output <- c(omega = omega, eta = eta, lambda = lambda)
  }

  if (from == 2 && to == 1) {
    chi <- delta^2
    psi <- gamma^2
    output <- c(chi = chi, psi = psi, lambda = lambda)
  }

  if (from == 2 && to == 3) {
    alpha <- gamma / delta
    beta <- gamma * delta
    output <- c(alpha = alpha, beta = beta, lambda = lambda)
  }

  if (from == 2 && to == 4) {
    omega <- delta * gamma
    eta <- delta / gamma
    output <- c(omega = omega, eta = eta, lambda = lambda)
  }

  if (from == 3 && to == 1) {
    chi <- beta / alpha
    psi <- alpha * beta
    output <- c(chi = chi, psi = psi, lambda = lambda)
  }

  if (from == 3 && to == 2) {
    delta <- sqrt(beta / alpha)
    gamma <- sqrt(alpha * beta)
    output <- c(delta = delta, gamma = gamma, lambda = lambda)
  }

  if (from == 3 && to == 4) {
    omega <- beta
    eta <- 1 / alpha
    output <- c(omega = omega, eta = eta, lambda = lambda)
  }

  if (from == 4 && to == 1) {
    chi <- omega * eta
    psi <- omega / eta
    output <- c(chi = chi, psi = psi, lambda = lambda)
  }

  if (from == 4 && to == 2) {
    delta <- sqrt(omega * eta)
    gamma <- sqrt(omega / eta)
    output <- c(delta = delta, gamma = gamma, lambda = lambda)
  }

  if (from == 4 && to == 3) {
    alpha <- 1 / eta
    beta <- omega
    output <- c(alpha = alpha, beta = beta, lambda = lambda)
  }

  if (from == to) {
    if (from == 1)
      output <- c(chi = chi, psi = psi, lambda = lambda)

    if (from == 2)
    output <- c(delta = delta, gamma = gamma, lambda = lambda)

    if (from == 3)
      output <- c(alpha = alpha, beta = beta, lambda = lambda)

    if (from == 4)
      output <- c(omega = omega, eta = eta, lambda = lambda)
  }

  if (noNames)
    names(output) <- NULL

  output
} ## End of gigChangePars
