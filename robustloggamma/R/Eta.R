### LOGGAMMA ETA=MEAN(EXP(Y)) PARAMETER

Exp.response <- function(mu, sigma, lambda, eps=0.0001, npoints=100000) {
  if (is.na(lambda)) return(NA)
  if (lambda > eps) {
    alpha <- 1/lambda^2
    ca <- lambda/sigma
    lMu <- mu - log(alpha)/ca + lgamma(alpha+1/ca)-lgamma(alpha)
    Mu <- exp(lMu)
  }
  if (abs(lambda) <= eps)
    Mu <- exp(mu + sigma^2/2)
  if (lambda < -eps) {
    ql <- qloggamma(p=ppoints(npoints),lambda=lambda)
    Mu <- mean(exp(mu+sigma*ql))
  }
  return(Mu)
}

### DERIVATA DI ETA RISPETTO A MU
Exp.response1 <- function(mu, sigma, lambda, eps=0.0001, npoints=100000) {
  Exp.response(mu, sigma, lambda, eps, npoints)
}

### DERIVATA DI ETA RISPETTO A SIGMA
Exp.response2 <- function(mu, sigma, lambda, eps=0.0001, npoints=100000) {
  if (is.na(lambda)) return(NA)
  Mu <- Exp.response(mu, sigma, lambda, eps, npoints)
  if (lambda > eps)
    Mu2 <- Mu*2*log(lambda)/lambda+Mu*digamma(lambda^(-2)+sigma/lambda)/lambda
  if (abs(lambda) <= eps)
    Mu2 <- Mu*sigma 
  if (lambda < -eps) {
    if (is.finite(Mu)) {
      myenv <- new.env()
      assign("x", sigma, envir = myenv)
      temp <- function(x) Exp.response(mu=mu, sigma=x, lambda=lambda, npoints=npoints)
      Mu2 <- drop(attr(numericDeriv(quote(temp(x)), "x", myenv), "gradient"))
    } else {
      Mu2 <- NaN
    }
  }
  return(Mu2)
}

### DERIVATA DI ETA RISPETTO A LAMBDA
Exp.response3 <- function(mu, sigma, lambda, eps=0.0001, npoints=100000) {
  if (is.na(lambda)) return(NA)
  Mu <- Exp.response(mu, sigma, lambda, eps, npoints)
  if (lambda > eps) {
    Mu3a <- 2*Mu*sigma*(1-log(lambda))/lambda^2
    Mu3b <- exp(mu + 2*log(lambda)*sigma/lambda)
    Mu3c <- gamma(lambda^(-2) + sigma/lambda)/gamma(lambda^(-2)) * (2*lambda^(-3)*digamma(lambda^(-2))-(2*lambda^(-3) + sigma*lambda^(-2))*digamma(lambda^(-2)+sigma/lambda))
    Mu3 <- Mu3a + Mu3b*Mu3c
  }
  if (abs(lambda) <= eps)
    Mu3 <- 0
  if (lambda < -eps) {
    if (is.finite(Mu)) {
      myenv <- new.env()
      assign("x", lambda, envir = myenv)
      temp <- function(x) Exp.response(mu=mu, sigma=sigma, lambda=x, npoints=npoints)
      Mu3 <- drop(attr(numericDeriv(quote(temp(x)), "x", myenv), "gradient"))       } else {
      Mu3 <- NaN
    }
  }
  return(Mu3)
}
