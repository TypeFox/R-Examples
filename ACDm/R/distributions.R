dburr <- function(x, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F){
  if(forceExpectation) theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
  
  retrunValue <- theta * kappa * x^(kappa - 1)
  retrunValue <- retrunValue / (1+ sig2 * theta * x^kappa)^(1/sig2 + 1)
  return(retrunValue)
}

pburr <- function(x, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F){
  if(forceExpectation) theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
  
  return(1 - (1 + sig2 * theta * x^kappa)^(-1/sig2))
}

qburr <- function(p, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F){
  if(forceExpectation) theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
  
  retrunValue <- (1-p)^(-sig2) - 1
  retrunValue <- retrunValue / (sig2 * theta)
  retrunValue <- retrunValue^(1/kappa)
  return(retrunValue)
}

burrExpectation <- function(theta = 1, kappa = 1.2, sig2 = .3){
  retrunValue <- theta^(-1/kappa) * gamma(1+1/kappa)*gamma(1/sig2-1/kappa)
  retrunValue <- retrunValue / (sig2^(1+1/kappa)*gamma(1/sig2+1))
  return(retrunValue)
}

rburr <- function(n = 1,theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F){
  if(forceExpectation) theta <- ((gamma(1+1/kappa)*gamma(1/sig2 - 1/kappa))/(sig2^(1+1/kappa)*gamma(1/sig2+1)))^(kappa)
  
  qburr(stats::runif(n), theta = theta, kappa = kappa, sig2 = sig2)
}

dgenf <- function(x, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F){
  if(min(kappa, eta, gamma, lambda) <= 0) stop("all parameters must be > 0")
  
  if(forceExpectation){
    if(eta - 1/gamma == 0) stop("The expectation is undefined for eta - 1/gamma == 0")
    
    lambda <- lgamma(kappa) + lgamma(eta) - lgamma(kappa + 1/gamma) - lgamma(eta - 1/gamma)
    lambda <- exp(lambda) * eta^(1/gamma)
  }
  
  returnValue <- (kappa * gamma -1) * log(x) + 
                 (-eta - kappa) * log(eta + (x/lambda)^gamma) + 
                 eta * log(eta) -
                 lbeta(kappa, eta) -
                 (kappa*gamma) * log(lambda)
  returnValue <- exp(returnValue) * gamma
  
  return(returnValue)  
}

pgenf <- function(q, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F){
  if(min(kappa, eta, gamma, lambda) <= 0) stop("all parameters must be > 0")
  
  if(forceExpectation){
    if(eta - 1/gamma == 0) stop("The expectation is undefined for eta - 1/gamma == 0")
    
    lambda <- lgamma(kappa) + lgamma(eta) - lgamma(kappa + 1/gamma) - lgamma(eta - 1/gamma)
    lambda <- exp(lambda) * eta^(1/gamma)
  }
  
  f <- function(x) dgenf(x = x, kappa = kappa, eta = eta, gamma = gamma, lambda = lambda)
  
  returnValue <- ifelse(q == 0, 0, stats::integrate(f, 0, q)$value)
  
  return(returnValue)  
}

genfHazard <- function(x, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F){
  
  if(min(kappa, eta, gamma, lambda) <= 0) stop("all parameters must be > 0")
  
  if(forceExpectation){
    if(eta - 1/gamma == 0) stop("The expectation is undefined for eta - 1/gamma == 0")
    
    lambda <- (lgamma(kappa) + lgamma(eta) 
              -(1/gamma)*log(eta)
              - lgamma(kappa + 1/gamma) - lgamma(eta - 1/gamma))
    lambda <- exp(lambda)
  }
    
  pdf <- dgenf(x = x, kappa = kappa, eta = eta, gamma = gamma, lambda = lambda)
  survivial <- 1 - pgenf(q = x, kappa = kappa, eta = eta, gamma = gamma, lambda)
  retrunValue <- pdf / survivial
  
  return(retrunValue)
  
}

dgengamma <- function(x, gamma = .3, kappa = 1.2, lambda = .3, forceExpectation = F){
  if(forceExpectation) lambda <- exp(lgamma(kappa) - lgamma(kappa + 1 / gamma))
  
  retrunValue <- ((kappa * gamma - 1) * log(x)
                  - (kappa * gamma) * log(lambda)
                  - lgamma(kappa)
                  + log(gamma)  
                  -(x / lambda)^gamma)
  
  retrunValue <- exp(retrunValue)
  
  return(retrunValue)
}

pgengamma <- function(x, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F){
  if(forceExpectation) lambda <- exp(lgamma(kappa) - lgamma(kappa + 1 / gamma))
    
  retrunValue <- stats::pgamma((x / lambda)^gamma, kappa)
  
  return(retrunValue)
}

qgengamma <- function(p, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F){
  if(forceExpectation) lambda <- exp(lgamma(kappa) - lgamma(kappa + 1 / gamma))
  
  retrunValue <- lambda * stats::qgamma(p, kappa)^(1/gamma)
  
  return(retrunValue)
}

rgengamma <- function(n = 1, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F){
  qgengamma(stats::runif(n), gamma = gamma, kappa = kappa, forceExpectation = forceExpectation)
}

gengammaHazard <- function(x, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F){
  if(forceExpectation) lambda <- exp(lgamma(kappa) - lgamma(kappa + 1 / gamma))
  
  pdf <- dgengamma(x, gamma = gamma, kappa = kappa, lambda = lambda)
  survivial <- 1 - pgengamma(x, gamma = gamma, kappa = kappa, lambda = lambda)
  retrunValue <- pdf / survivial
  
  return(retrunValue)
}

dqweibull <- function(x, a = .8, qdist = 1.2, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
    
  returnValue <- (2 - qdist) * a / b^a * x^(a-1) * (1-(1-qdist)*(x/b)^a)^(1/(1-qdist)) 
  
  return(returnValue)
}

pqweibull <- function(q, a = .8, qdist = 1.2, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
    
  returnValue <- 1 - (1-(1-qdist)*(q/b)^a)^((2-qdist)/(1-qdist))
  
  return(returnValue)
}

qqweibull <- function(p, a = .8, qdist = 1.2, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 6, distPara = c(a, q))
  
  returnValue <- ((1 - (1 - p)^((1 - qdist) / (2 - qdist))) / (1 - qdist))^(1 / a) * b  
  return(returnValue)
}

rqweibull <- function(n = 1, a = .8, qdist = 1.2, b = 1, forceExpectation = F){
  
  qqweibull(stats::runif(n), a = a, qdist = qdist, b = b, forceExpectation = forceExpectation)
  
}

qweibullExpectation <- function(a = .8, qdist = 1.2, b = 1){
  
  if((1/(qdist-1)-1)*a <= 1) stop("expectation does not exist for the given parameters")
  
  returnValue <- lgamma(1/a) + lgamma(1/(qdist-1) - 1/a - 1) - lgamma(1/(qdist-1))
  
  returnValue <- exp(returnValue) * (2 - qdist) / ((qdist - 1)^((a+1)/a) * a)
  
  returnValue <- returnValue * b
  
  return(returnValue)
}

qweibullHazard <- function(x, a = .8, qdist = 1.2, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
  
  returnValue <- (2-qdist) * b^(-a) * x^(a-1) * a
  returnValue <- returnValue / (1 - (1-qdist) * (x/b)^a)
  
  return(returnValue)
}

dmixqwe <- function(x, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F){
  if(forceExpectation) b <- (1 - (1 - pdist) * lambda) / pdist * .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
  
  returnValue <- pdist * (2 - qdist) * a / b^a * x^(a-1) * (1-(1-qdist)*(x/b)^a)^(1/(1-qdist)) 
  returnValue <- returnValue + (1 - pdist) * 1/lambda * exp(-x/lambda)
  
  return(returnValue)
}

pmixqwe <- function(q, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F){
  if(forceExpectation) b <- (1 - (1 - pdist) * lambda) / pdist * .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
  
  returnValue <- pdist * pqweibull(q, a = a, q = qdist, b = b) + (1 - pdist) * (1 - exp(-q/lambda))
  
  return(returnValue)
}

mixqweHazard <- function(x, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F){
  if(forceExpectation) b <- (1 - (1 - pdist) * lambda) / pdist * .returnFixedMeanPara(distCode = 6, distPara = c(a, qdist))
  
  pdf <- dmixqwe(x, pdist = pdist, a = a, qdist = qdist, lambda = lambda, b = b)
  survivial <- 1 - pmixqwe(x, pdist = pdist, a = a, qdist = qdist, lambda = lambda, b = b)
  returnValue <- pdf / survivial
  
  return(returnValue)
}

dmixqww <- function(x, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 8, distPara = c(pdist, a, qdist, theta, gamma))
  
  returnValue <- pdist * (2 - qdist) * a / b^a * x^(a-1) * (1-(1-qdist)*(x/b)^a)^(1/(1-qdist)) 
  returnValue <- returnValue + (1 - pdist) * theta * gamma * x^(gamma-1) * exp(-theta * x^gamma)
  
  return(returnValue)
}

pmixqww <- function(q, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 8, distPara = c(pdist, a, qdist, theta, gamma))
    
  returnValue <- pdist * pqweibull(q, a = a, q = qdist, b = b) + (1 - pdist) * (1 - exp(-theta * q^gamma))
  
  return(returnValue)
}

mixqwwHazard <- function(x, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1, forceExpectation = F){
  if(forceExpectation) b <- .returnFixedMeanPara(distCode = 8, distPara = c(pdist, a, qdist, theta, gamma))
  
  pdf <- dmixqww(x, pdist = pdist, a = a, qdist = qdist, theta = theta, gamma = gamma, b = b, forceExpectation = forceExpectation)
  survivial <- 1 - pmixqww(x, pdist = pdist, a = a, qdist = qdist, theta = theta, gamma = gamma, b = b, forceExpectation = forceExpectation)
  returnValue <- pdf / survivial
  
  return(returnValue)
}


dmixinvgauss <- function(x, theta = .2, lambda = .1, gamma = .05, forceExpectation = F){
  phi <- 1
  if(forceExpectation) phi <- theta * (1 + theta^2 / lambda / (1 + gamma))
  
  returnValue <- (gamma + phi * x) / (gamma + theta) * sqrt(lambda/(2 * pi * x^3 * phi)) * exp(-lambda * (phi * x - theta)^2/(2 * phi * x * theta^2))
  
  return(returnValue)
}


pmixinvgauss <- function(q, theta = .2, lambda = .1, gamma = .05, forceExpectation = F){
  phi <- 1
  if(forceExpectation) phi <- theta * (1 + theta^2 / lambda / (1 + gamma))
  
  t1 <- q / theta - 1
  t2 <- -q / theta - 1
  
  returnValue <- stats::pnorm(t1*sqrt(lambda / (q * phi))) + (gamma - theta) / (theta + gamma) * stats::pnorm(t2*sqrt(lambda / (q * phi))) * exp(2 * lambda / theta)
  
  return(returnValue)
}

mixinvgaussHazard <- function(x, theta = .2, lambda = .1, gamma = .05, forceExpectation = F){
  
  pdf <- dmixinvgauss(x, theta = theta, lambda = lambda, gamma = gamma, forceExpectation = forceExpectation)
  survivial <- 1 - pmixinvgauss(x, theta = theta, lambda = lambda, gamma = gamma, forceExpectation = forceExpectation)
  returnValue <- pdf / survivial
  
  return(returnValue)
}
# 
# dbirnbaumsaunders <- function(x, kappa = 1, sigma = 1){
#   
#   returnValue <- 1/(2 * kappa * sigma * sqrt(2 * pi)) * ((sigma/x)^(1/2) + (sigma/x)^(3/2)) * exp(-1/(2 * kappa^2) * (x/sigma + sigma/x -2))
#   
#   return(returnValue)
# }
# 
# pbirnbaumsaunders <- function(x, kappa = 1, sigma = 1){
#   
#   returnValue <- stats::pnorm(1/kappa * ((x/sigma)^(1/2) - (sigma/x)^(1/2)))
#   
#   return(returnValue)
# }
# 
# birnbaumsaundersHazard <- function(x, kappa = 1, sigma = 1){
#   
#   pdf <- dbirnbaumsaunders(x = x, kappa = kappa, sigma = sigma)
#   survivial <- 1 - pbirnbaumsaunders(x, kappa, sigma = sigma)
#   returnValue <- pdf / survivial
#   
#   return(returnValue)
# }