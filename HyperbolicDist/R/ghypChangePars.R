### Change parameterizations of the Generalized Hyperbolic Distribution
ghypChangePars <- function (from, to, Theta, noNames = FALSE) 
{
  Theta <- as.numeric(Theta)
  if(length(Theta)==4) {
    Theta <- c(1,Theta)
  }else{
    if (length(Theta) != 5) {
      stop("parameter vector must contain 4 or 5 values")
    }
  }
  if ((from != 1) & (from != 2) & (from != 3) & (from != 4)) {
    stop("the argument 'from' must be either 1, 2, 3 or 4")
  }
  if ((to != 1) & (to != 2) & (to != 3) & (to != 4)) {
    stop("the argument 'to' must be either 1, 2, 3 or 4")
  }
  lambda <- Theta[1]
  delta <- Theta[4]
  if (delta <= 0) 
    stop("delta must be greater than zero")
  mu <- Theta[5]
  if (from == 1) {
     alpha <- Theta[2]
     beta <- Theta[3]
    if (alpha <= 0) 
      stop("alpha must be greater than zero")
    if (abs(beta) >= alpha) 
      stop("absolute value of beta must be less than alpha")
  }
  if (from == 2) {
    zeta <- Theta[2]
    rho <- Theta[3]
    if (zeta <= 0) 
      stop("zeta must be greater than zero")
  }
  if (from == 3) {
    xi <- Theta[2]
    chi <- Theta[3]
    if ((xi <= 0) || (xi >1)) 
      stop("xi must be between zero and one")
    if (abs(chi) > xi) 
      stop("absolute value of chi must be less than xi")
  }
  if (from == 4) {
    alphaBar <- Theta[2]
    betaBar <- Theta[3]
    if (alphaBar <= 0) 
      stop("alpha bar must be greater than zero")
    if (abs(betaBar) >= alphaBar) 
      stop("absolute value of beta bar must be less than alpha bar")
  }
  if (from == 1 && to == 2) {
    zeta <- delta*sqrt(alpha^2 - beta^2)
    rho <- beta/alpha
    output <- c(lambda = lambda, zeta = zeta, rho = rho,
                delta = delta, mu = mu)
  }
  if (from == 1 && to == 3) {
    xi <- 1/sqrt(1 + delta*sqrt(alpha^2 - beta^2))
    chi <- beta/(alpha*sqrt(1 + delta*sqrt(alpha^2 - beta^2)))
    output <- c(lambda = lambda, xi = xi, chi = chi,
                delta = delta, mu = mu)
  }
  if (from == 1 && to == 4) {
    alphaBar <- delta*alpha
    betaBar <- delta*beta
    output <- c(lambda = lambda, alphaBar = alphaBar, betaBar = betaBar,
                delta = delta, mu = mu)
  }
  if (from == 2 && to == 1) {
    alpha <- zeta/(delta*sqrt(1 - rho^2))
    beta <- rho*alpha
    output <- c(lambda = lambda, alpha = alpha, beta = beta,
                delta = delta, mu = mu)
  }
  if (from == 2 && to == 3) {
    xi <- 1/sqrt(1 + zeta)
    chi <- xi*rho
    output <- c(lambda = lambda, xi = xi, chi = chi,
                delta = delta, mu = mu)
  }
  if (from == 2 && to == 4) {
    alphaBar <- zeta/sqrt(1 - rho^2)
    betaBar <- rho*alphaBar
    output <- c(lambda = lambda, alphaBar = alphaBar, betaBar = betaBar,
                delta = delta, mu = mu)
  }
  if (from == 3 && to == 1) {
    alpha <- (1 - xi^2)/(delta*xi*sqrt(xi^2 - chi^2))
    beta <- alpha*chi/xi
    output <- c(lambda = lambda, alpha = alpha, beta = beta,
                delta = delta, mu = mu)
  }
  if (from == 3 && to == 2) {
    zeta <- (1/xi^2) - 1
    rho <- chi/xi
    output <- c(lambda = lambda, zeta = zeta, rho = rho,
                delta = delta, mu = mu)
  }
  if (from == 3 && to == 4) {
    alphaBar <- (1 - xi^2)/(xi*sqrt(xi^2 - chi^2))
    betaBar <- alphaBar*chi/xi
    output <- c(lambda = lambda, alphaBar = alphaBar, betaBar = betaBar,
                delta = delta, mu = mu)
  }
  if (from == 4 && to == 1) {
    alpha <- alphaBar/delta
    beta <- betaBar/delta
    output <- c(lambda = lambda, alpha = alpha, beta = beta,
                delta = delta, mu = mu)
  }
  if (from == 4 && to == 2) {
    zeta <- sqrt(alphaBar^2 - betaBar^2)
    rho <- betaBar/alphaBar
    output <- c(lambda = lambda, zeta = zeta, rho = rho,
                delta = delta, mu = mu)
  }
  if (from == 4 && to == 3) {
    xi <- 1/sqrt(1 + sqrt(alphaBar^2 - betaBar^2))
    chi <- betaBar/(alphaBar*sqrt(1 + sqrt(alphaBar^2 - betaBar^2)))
    output <- c(lambda = lambda, xi = xi, chi = chi,
                delta = delta, mu = mu)
  }
  if (from == to) {
    if (from == 1) 
      output <- c(lambda = lambda, alpha = alpha, beta = beta,
                  delta = delta, mu = mu)
    if (from == 2) 
      output <- c(lambda = lambda, zeta = zeta, rho = rho,
                  delta = delta, mu = mu)
    if (from == 3) 
      output <- c(lambda = lambda, xi = xi, chi = chi,
                  delta = delta, mu = mu)
    if (from == 4) 
      output <- c(lambda = lambda, alphaBar = alphaBar, betaBar = betaBar,
                  delta = delta, mu = mu)
  }
  if (noNames == TRUE) 
      names(output) <- NULL
  output
} ## End of hyperbChangePars()
