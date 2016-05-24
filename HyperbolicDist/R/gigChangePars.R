### Change parameterizations of the generalized inverse Gaussian distribution
gigChangePars <- function (from, to, Theta, noNames = FALSE) 
{
  if (length(Theta) != 3) {
    stop("parameter vector must contain 3 values")
  }
  if ((from != 1) & (from != 2) & (from != 3) & (from != 4)) {
    stop("the argument 'from' must be either 1, 2, 3 or 4")
  }
  if ((to != 1) & (to != 2) & (to != 3) & (to != 4)) {
    stop("the argument 'to' must be either 1, 2, 3 or 4")
  }
  lambda <- Theta[1]
  if (from == 1) {
    chi <- Theta[2]
    psi <- Theta[3]
    if (chi <= 0) 
      stop("chi must be greater than zero")
    if (psi <= 0) 
      stop("psi must be greater than zero")
  }
  if (from == 2) {
    delta <- Theta[2]
    gamma <- Theta[3]
    if (delta <= 0) 
      stop("delta must be greater than zero")
    if (gamma <= 0) 
      stop("gamma must be greater than zero")
  }
  if (from == 3) {
    alpha <- Theta[2]
    beta <- Theta[3]
    if (alpha <= 0) 
      stop("alpha must be greater than zero")
     if (beta <= 0) 
       stop("beta must be greater than zero")
  }
  if (from == 4) {
    omega <- Theta[2]
    eta <- Theta[3]
    if (omega <= 0) 
      stop("omega must be greater than zero")
    if (eta <= 0) 
      stop("eta must be greater than zero")
  }
  if (from == 1 && to == 2) {
    delta <- sqrt(chi)
    gamma <- sqrt(psi)
    output <- c(lambda = lambda, delta = delta, gamma = gamma)
  }
  if (from == 1 && to == 3) {
    alpha <- sqrt(psi/chi)
    beta <- sqrt(chi*psi)
    output <- c(lambda = lambda, alpha = alpha, beta = beta)
  }
  if (from == 1 && to == 4) {
    omega <- sqrt(chi*psi)
    eta <- sqrt(chi/psi)
    output <- c(lambda = lambda, omega = omega, eta = eta)
  }
  if (from == 2 && to == 1) {
    chi <- delta^2
    psi <- gamma^2
    output <- c(lambda = lambda, chi = chi, psi = psi)
  }
  if (from == 2 && to == 3) {
    alpha <- gamma/delta
    beta <- gamma*delta
    output <- c(lambda = lambda, alpha = alpha, beta = beta)
  }
  if (from == 2 && to == 4) {
    omega <- delta*gamma
    eta <- delta/gamma
    output <- c(lambda = lambda, omega = omega, eta = eta)
  }
  if (from == 3 && to == 1) {
    chi <- beta/alpha
    psi <- alpha*beta
    output <- c(lambda = lambda, chi = chi, psi = psi)
  }
  if (from == 3 && to == 2) {
    delta <- sqrt(beta/alpha)
    gamma <- sqrt(alpha*beta)
    output <- c(lambda = lambda, delta = delta, gamma = gamma)
  }
  if (from == 3 && to == 4) {
    omega <- beta
    eta <- 1/alpha
    output <- c(lambda = lambda, omega = omega, eta = eta)
  }
  if (from == 4 && to == 1) {
    chi <- omega*eta
    psi <- omega/eta
    output <- c(lambda = lambda, chi = chi, psi = psi)
  }
  if (from == 4 && to == 2) {
    delta <- sqrt(omega*eta)
    gamma <- sqrt(omega/eta)
    output <- c(lambda = lambda, delta = delta, gamma = gamma)
  }
  if (from == 4 && to == 3) {
    alpha <- 1/eta
    beta <- omega
    output <- c(lambda = lambda, alpha = alpha, beta = beta)
  }
  if (from == to) {
    if (from == 1) 
      output <- c(lambda = lambda, chi = chi, psi = psi)
    if (from == 2) 
    output <- c(lambda = lambda, delta = delta, gamma = gamma)
    if (from == 3) 
      output <- c(lambda = lambda, alpha = alpha, beta = beta)
    if (from == 4) 
      output <- c(lambda = lambda, omega = omega, eta = eta)
  }
  if (noNames == TRUE) 
    names(output) <- NULL
  output
} ## End of gigChangePars
