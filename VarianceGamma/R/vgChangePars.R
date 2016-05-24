vgChangePars <- function (from, to, param, noNames = FALSE) {
  if (length(param)!=4) {
    stop("parameter vector must contain 4 values")
  }
  if ((from != 1) & (from != 2) & (from != 3) & (from != 4)) {
    stop("the argument 'from' must be either 1, 2, 3 or 4")
  }
  if ((to != 1) & (to != 2) & (to != 3) & (to != 4)) {
    stop("the argument 'to' must be either 1, 2, 3 or 4")
  }

  if (from == 1) {
    vgC <- param[1]
    sigma1 <- param[2]
    vgTheta1 <- param[3]
    nu <- param[4]
    if (nu <= 0)
      stop ("nu must be greater than zero")
    if (sigma1 <= 0)
      stop ("sigma must be greater than zero")
  }

  if (from == 2) {
    vgTheta2 <- param[1]
    sigma2 <- param[2]
    mu <- param[3]
    tau2 <- param[4]
    if (tau2 <= 0)
      stop ("tau must be greater than zero")
    if (sigma2 <= 0)
      stop ("sigma must be greater than zero")
  }

  if (from == 3) {
    vgTheta3 <- param[1]
    sigma3 <- param[2]
    kappa <- param[3]
    tau3 <- param[4]
    if (tau3 <= 0)
      stop ("tau must be greater than zero")
    if (sigma3 <= 0)
      stop ("sigma must be greater than zero")
  }

  if(from == 4) {
    lambda <- param[1]
    alpha <- param[2]
    beta <- param[3]
    mu <- param[4]
    if (lambda <= 0)
      stop ("lambda must be greater than zero")
    if (alpha <= 0) {
      stop("alpha must be greater than zero")
    }
    if (abs(beta) >= alpha) {
      stop("absolute value of beta must be less than alpha")
    }

  }

  if (from == 1 && to == 2) {
    mu <- vgTheta1*nu
    tau2 <- 1/nu
    sigma2 <- sigma1*sqrt(nu)
    vgTheta2 <- vgC
    output = c(theta = vgTheta2,sigma = sigma2,mu = mu,tau = tau2)
  }

  if (from == 1 && to == 3) {
    kappa <- (sqrt(2*(sigma1*sqrt(nu))^2 + (vgTheta1*nu)^2) - vgTheta1*nu)/
      (sigma1*sqrt(nu)*sqrt(2))
    tau3 <- 1/nu
    sigma3 <- sigma1*sqrt(nu)
    vgTheta3 <- vgC
    output = c(theta = vgTheta3,sigma = sigma3,kappa = kappa,tau = tau3)
  }

  if (from == 1 && to == 4) {
    lambda <- 1/nu
    beta <- vgTheta1/sigma1^2
    alpha <- sqrt(2/(nu*sigma1^2) + vgTheta1^2/sigma1^4)
    mu <- vgC
    output = c(lambda = lambda,alpha = alpha,beta = beta,mu = mu)
  }

  if (from == 2 && to == 1) {
    vgTheta1 <- mu*tau2
    nu <- 1/tau2
    sigma1 <- sqrt(tau2)*sigma2
    vgC <- vgTheta2
    output = c(c = vgC,sigma = sigma1,theta = vgTheta1,nu = nu)
  }

  if (from == 2 && to == 3) {
    kappa <- (sqrt(2*(sigma2)^2 + mu^2) - mu)/(sigma2*sqrt(2))
    tau3 <- tau2
    sigma3 <- sigma2
    vgTheta3 <-vgTheta2
    output = c(theta = vgTheta3,sigma = sigma3,kappa = kappa,tau = tau3)
  }

  if (from == 2 && to == 4) {
    vgTheta1 <- mu*tau2
    nu <- 1/tau2
    sigma1 <- sqrt(tau2)*sigma2
    vgC <- vgTheta2
    lambda <- 1/nu
    beta <- vgTheta1/sigma1^2
    alpha <- sqrt(2/(nu*sigma1^2) + vgTheta1^2/sigma1^4)
    mu <- vgC
    output = c(lambda = lambda,alpha = alpha,beta = beta,mu = mu)
  }

  if (from == 3 && to == 1) {
    vgTheta1 <- tau3*sigma3*(1/kappa - kappa)/sqrt(2)
    nu <- 1/tau3
    sigma1 <- sqrt(tau3)*sigma3
    vgC <- vgTheta3
    output = c(c = vgC,sigma = sigma1,theta = vgTheta1,nu = nu)
  }

  if (from == 3 && to == 2) {
    mu <- sigma3*(1/kappa - kappa)/sqrt(2)
    tau2 <- tau3
    sigma2 <- sigma3
    vgTheta2 <- vgTheta3
    output = c(theta= vgTheta2,sigma = sigma2,mu = mu,tau = tau2)
  }

  if (from == 3 && to == 4) {
    vgTheta1 <- tau3*sigma3*(1/kappa - kappa)/sqrt(2)
    nu <- 1/tau3
    sigma1 <- sqrt(tau3)*sigma3
    vgC <- vgTheta3
    lambda <- 1/nu
    beta <- vgTheta1/sigma1^2
    alpha <- sqrt(2/(nu*sigma1^2) + vgTheta1^2/sigma1^4)
    mu <- vgC
    output = c(lambda = lambda,alpha = alpha,beta = beta,mu = mu)
  }

  if (from == 4 && to == 1) {
    sigma1 <- sqrt(2*lambda/(alpha^2 - beta^2))
    nu <- 1/lambda
    vgTheta1 <- beta*sigma1^2
    vgC <- mu
    output = c(c = vgC,sigma = sigma1,theta = vgTheta1,nu = nu)
  }

  if (from == 4 && to == 2) {
    sigma1 <- sqrt(2*lambda/(alpha^2 - beta^2))
    nu <- 1/lambda
    vgTheta1 <- beta*sigma1^2
    vgC <- mu
    mu <- vgTheta1*nu
    tau2 <- 1/nu
    sigma2 <- sigma1*sqrt(nu)
    vgTheta2 <- vgC
    output = c(theta = vgTheta2,sigma = sigma2,mu = mu,tau = tau2)
  }

  if (from == 4 && to == 3) {
    sigma1 <- sqrt(2*lambda/(alpha^2 - beta^2))
    nu <- 1/lambda
    vgTheta1 <- beta*sigma1^2
    vgC <- mu
    kappa <- (sqrt(2*(sigma1*sqrt(nu))^2 + (vgTheta1*nu)^2) - vgTheta1*nu)/
      (sigma1*sqrt(nu)*sqrt(2))
    tau3 <- 1/nu
    sigma3 <- sigma1*sqrt(nu)
    vgTheta3 <- vgC
    output = c(theta = vgTheta3,sigma = sigma3,kappa = kappa,tau = tau3)
  }

  if (from == to) {
    if (from == 1)
      output = c(c = vgC,sigma = sigma1,theta = vgTheta1,nu = nu)
    if (from == 2)
      output = c(theta = vgTheta2,sigma = sigma2,mu = mu,tau = tau2)
    if (from == 3)
      output = c(theta = vgTheta3,sigma = sigma3,kappa = kappa,tau = tau3)
    if (from == 4)
      output == c(lambda = lambda,alpha = alpha,beta = beta,mu = mu)
  }

  if (noNames == TRUE)
    names(output) <- NULL
  output
}

