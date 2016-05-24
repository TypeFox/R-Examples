### < ====================================================================== >
setGeneric("pstate",
           function(lower, upper, time = 1, s0, ...)
           standardGeneric("pstate"))


### < ---------------------------------------------------------------------- >
pstate.default <- function(lower, upper, time = 1, s0 = 50, delta0 = 0,
                           mu = 0.1, sigmaS = 0.3, kappa = 1, alpha = 0,
                           sigmaE = 0.5, rho = 0.75, ...)
{
  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, time = time)

  lower[1] <- log(lower[1])
  upper[1] <- log(upper[1])

  mean <- .mu.state.schwartz2f(x0 = log(s0), delta0 = delta0,
                               mu = mu, sigmaS = sigmaS,
                               kappa = kappa, alpha = alpha,
                               sigmaE = sigmaE, rho = rho,
                               time = time)

  sigma <- .sigma.state.schwartz2f(sigmaS = sigmaS, kappa = kappa,
                                   sigmaE = sigmaE, rho = rho, time = time)

  return(pmvnorm(lower = lower, upper = upper, mean = mean,
                 sigma = sigma, ...))
}

setMethod("pstate", signature(lower = "ANY", upper = "ANY", time = "ANY",
                              s0 = "numeric"),
          pstate.default)

### < ---------------------------------------------------------------------- >
pstate.schwartz2f <- function(lower, upper, time = 1, s0, ...)
{
  tmp.coef <- coef(s0)
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(pstate.default(lower, upper, time = time, s0 = tmp.coef$s0,
                        delta0 = delta0, mu = mu, sigmaS = sigmaS,
                        kappa = kappa, alpha = alpha, sigmaE = sigmaE,
                        rho = rho, ...))
}

setMethod("pstate", signature(lower = "ANY", upper = "ANY", time = "ANY",
                              s0 = "schwartz2f"),
          pstate.schwartz2f)
### < ---------------------------------------------------------------------- >


### < ====================================================================== >
setGeneric("dstate",
           function(x, time = 1, s0, ...)
           standardGeneric("dstate"))

dstate.default <- function(x, time = 1, s0 = 50, delta0 = 0, mu = 0.1,
                           sigmaS = 0.3, kappa = 1, alpha = 0,
                           sigmaE = 0.5, rho = 0.75, ...)
{

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, time = time)

 if(is.vector(x)){
    x <- rbind(x)
  }

  x <- as.matrix(x) # in case x is a data.frame

  x[,1] <- log(x[,1])

  mean <- .mu.state.schwartz2f(x0 = log(s0), delta0 = delta0,
                               mu = mu, sigmaS = sigmaS,
                               kappa = kappa, alpha = alpha,
                               sigmaE = sigmaE, rho = rho,
                               time = time)

  sigma <- .sigma.state.schwartz2f(sigmaS = sigmaS, kappa = kappa,
                                   sigmaE = sigmaE, rho = rho, time = time)

  dens <- dmvnorm(x, mean = mean, sigma = sigma) / exp(x[,1], ...)

  return(dens)

}


setMethod("dstate", signature(x = "ANY", time = "ANY", s0 = "numeric"),
          dstate.default)
### < ---------------------------------------------------------------------- >

dstate.schwartz2f <- function(x, time = 1, s0, ...)
{
  tmp.coef <- coef(s0)
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(dstate.default(x, time = time, s0 = tmp.coef$s0, delta0 = delta0,
                        mu = mu, sigmaS = sigmaS, kappa = kappa,
                        alpha = alpha, sigmaE = sigmaE, rho = rho, ...))
}

setMethod("dstate", signature(x = "ANY", time = "ANY",
                              s0 = "schwartz2f"),
          dstate.schwartz2f)
### < ---------------------------------------------------------------------- >

### < ====================================================================== >
setGeneric("qstate",
           function(p, time = 1, s0, ...)
           standardGeneric("qstate"))


### < ---------------------------------------------------------------------- >
qstate.default <- function(p, time = 1, s0 = 50, delta0 = 0, mu = 0.1,
                           sigmaS = 0.3, kappa = 1, alpha = 0,
                           sigmaE = 0.5, rho = 0.75,
                           tail = "lower.tail", ...)
{
  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, time = time)

  mean <- .mu.state.schwartz2f(x0 = log(s0), delta0 = delta0,
                               mu = mu, sigmaS = sigmaS,
                               kappa = kappa, alpha = alpha,
                               sigmaE = sigmaE, rho = rho,
                               time = time)

  sigma <- .sigma.state.schwartz2f(sigmaS = sigmaS, kappa = kappa,
                                   sigmaE = sigmaE, rho = rho, time = time)

  quant <- qmvnorm(p = p, tail = tail, mean = mean, sigma = sigma, ...)
  quant$quantile <- c(exp(quant$quantile), quant$quantile)

  return(quant)
}

setMethod("qstate", signature(p = "ANY", time = "ANY", s0 = "numeric"),
          qstate.default)
### < ---------------------------------------------------------------------- >

qstate.schwartz2f <- function(p, time = 1, s0, tail = "lower.tail", ...)
{
  tmp.coef <- coef(s0)
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(qstate.default(p, time = time, s0 = tmp.coef$s0, delta0 = delta0,
                        mu = mu, sigmaS = sigmaS,
                        kappa = kappa, alpha = alpha,
                        sigmaE = sigmaE, rho = rho,
                        tail = tail, ...))
}

setMethod("qstate", signature(p = "ANY", time = "ANY", s0 = "schwartz2f"),
          qstate.schwartz2f)
### < ---------------------------------------------------------------------- >

### < ====================================================================== >
setGeneric("rstate",
           function(n, time = 1, s0, ...)
           standardGeneric("rstate"))

### < ---------------------------------------------------------------------- >
rstate.default <- function(n, time = 1, s0 = 50, delta0 = 0,
                           mu = 0.1, sigmaS = 0.3,
                           kappa = 1, alpha = 0,
                           sigmaE = 0.5, rho = 0.75,
                           method = "chol")
{
  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, time = time)
  
  if(n != round(n)){
    stop("'n' must be an integer number!")
  }

  mean <- .mu.state.schwartz2f(x0 = log(s0), delta0 = delta0,
                               mu = mu, sigmaS = sigmaS,
                               kappa = kappa, alpha = alpha,
                               sigmaE = sigmaE, rho = rho,
                               time = time)

  sigma <- .sigma.state.schwartz2f(sigmaS = sigmaS, kappa = kappa,
                                   sigmaE = sigmaE, rho = rho, time = time)

  rand <- rmvnorm(n = n, mean = mean, sigma = sigma, method = method)

##  browser()
  rand[,1] <- exp(rand[,1])
##  browser()
  colnames(rand) <- c("S", "delta")
  return(rand)
}

setMethod("rstate", signature(n = "ANY", time = "ANY", s0 = "numeric"),
          rstate.default)
### < ---------------------------------------------------------------------- >

rstate.schwartz2f <- function(n, time = 1, s0, method = "chol")
{
  tmp.coef <- coef(s0)
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(rstate.default(n, time = time, s0 = tmp.coef$s0, delta0 = delta0,
                        mu = mu, sigmaS = sigmaS,
                        kappa = kappa, alpha = alpha,
                        sigmaE = sigmaE, rho = rho,
                        method = method))
}

setMethod("rstate", signature(n = "ANY", time = "ANY", s0 = "schwartz2f"),
          rstate.schwartz2f)
### < ---------------------------------------------------------------------- >


### < ====================================================================== >
setGeneric("simstate",
           function(n, time = 1, s0, ...)
           standardGeneric("simstate"))

### < ---------------------------------------------------------------------- >
simstate.default <- function(n, time = 1, s0 = 50, delta0 = 0,
                             mu = 0.1, sigmaS = 0.3,
                             kappa = 1, alpha = 0,
                             sigmaE = 0.5, rho = 0.75,
                             method = "chol")
{
  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, time = time)

  if(n <= 1){
    stop("Use 'rstate' if n <= 1!")
  }

  deltat <- time/n

  sigma <- .sigma.state.schwartz2f(sigmaS = sigmaS, kappa = kappa,
                                   sigmaE = sigmaE, rho = rho, time = deltat)

  traj <- matrix(NA, ncol = 2, nrow = n,
                 dimnames = list(1:n, c("S", "delta")))
  traj[1, ] <- c(log(s0), delta0)

  increments <- rmvnorm(n = n - 1, mean = c(0, 0),
                        sigma = sigma, method = method)

  for(i in 2:n){
    drift <- .mu.state.schwartz2f(x0 = traj[i - 1, 1],
                                  delta0 = traj[i - 1, 2],
                                  mu = mu, sigmaS = sigmaS,
                                  kappa = kappa, alpha = alpha,
                                  sigmaE  = sigmaE, rho = rho,
                                  time = deltat, as.mat = FALSE)
    traj[i, ] <- drift + increments[i - 1, ]
  }

  traj[, 1] <- exp(traj[, 1])

  return(traj)
}

setMethod("simstate", signature(n = "ANY", time = "ANY", s0 = "numeric"),
          simstate.default)
### < ---------------------------------------------------------------------- >

simstate.schwartz2f <- function(n, time = 1, s0, method = "chol")
{
  tmp.coef <- coef(s0)
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(simstate.default(n, time = time, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa,
                          alpha = alpha, sigmaE = sigmaE, rho = rho,
                          method = method))

}

setMethod("simstate", signature(n = "ANY", time = "ANY", s0 = "schwartz2f"),
          simstate.schwartz2f)
### < ---------------------------------------------------------------------- >

### < ====================================================================== >
setGeneric("pfutures",
           function(q, time = 0.1, ttm = 1, s0, ...)
           standardGeneric("pfutures"))

### < ---------------------------------------------------------------------- >
pfutures.default <- function(q, time = 0.1, ttm = 1, s0 = 50, delta0 = 0, mu = 0.1,
                             sigmaS = 0.3, kappa = 1, alpha = 0,
                             sigmaE = 0.5, rho = 0.75,
                             r = 0.05, lambda = 0, alphaT = NULL,
                             measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)
  
  stopifnot(time <= ttm)
  if((missing(lambda) | missing(alpha)) & missing(alphaT)){
    warning("Both 'alphaT' and ('lambda' or 'alpha') are missing!\n",
            "The mean-level of the convenience yield is set to zero.")
    alphaT <- 0
  }else if(missing(alphaT)){
    alphaT <- alpha - lambda / kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda, alphaT = alphaT,
                 r = r, time = time, ttm = ttm)
  
  mu.fut <- .mu.fut.schwartz2f(x0 = log(s0), delta0 = delta0, mu = mu,
                               sigmaS = sigmaS, kappa = kappa,
                               sigmaE = sigmaE, rho = rho, alpha = alpha,
                               alphaT = alphaT, r = r, time = time, ttm = ttm,
                               measure = measure)

  sigma.fut <- .sigma.fut.schwartz2f(sigmaS = sigmaS, kappa = kappa, sigmaE = sigmaE,
                                     rho = rho, time = time, ttm = ttm)

  return(pnorm(log(q), mean = mu.fut, sd = sqrt(sigma.fut), ...))
}

setMethod("pfutures", signature(q = "ANY", time ="ANY", ttm = "ANY", s0 = "numeric"),
          pfutures.default)
### < ---------------------------------------------------------------------- >

pfutures.schwartz2f <- function(q, time = 0.1, ttm = 1, s0, r = 0.05, lambda = 0,
                                alphaT = NULL,
                                measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)

  tmp.coef <- coef(s0)

  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    alphaT <- tmp.coef$alpha
  }else if(missing(alphaT)){
    alphaT <- tmp.coef$alpha - lambda / tmp.coef$kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(pfutures.default(q, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho, r = r,
                          alphaT = alphaT, measure = measure, ...))
}

setMethod("pfutures", signature(q = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f"),
          pfutures.schwartz2f)
### < ---------------------------------------------------------------------- >

pfutures.schwartz2f.fit <- function(q, time = 0.1, ttm = 1, s0,
                                    measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)

  tmp.coef <- coef(s0)

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  alpha <- tmp.coef$alpha
  kappa <- tmp.coef$kappa
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho
  r <- tmp.coef$r
  alphaT <- tmp.coef$alphaT

  return(pfutures.default(q, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho,
                          r = r, alphaT = alphaT, measure = measure, ...))
}

setMethod("pfutures", signature(q = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f.fit"),
          pfutures.schwartz2f.fit)
### < ---------------------------------------------------------------------- >

### < ====================================================================== >
setGeneric("dfutures",
           function(x, time = 0.1, ttm = 1, s0, ...)
           standardGeneric("dfutures"))

dfutures.default <- function(x, time = 0.1, ttm = 1, s0 = 50, delta0 = 0, mu = 0.1,
                             sigmaS = 0.3, kappa = 1, alpha = 0,
                             sigmaE = 0.5, rho = 0.75, r = 0.05,
                             lambda = 0, alphaT = NULL,
                             measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)

  stopifnot(time <= ttm)
  if((missing(lambda) | missing(alpha)) & missing(alphaT)){
    warning("Both 'alphaT' and ('lambda' or 'alpha') are missing!\n",
            "The mean-level of the convenience yield is set to zero.")
    alphaT <- 0
  }else if(missing(alphaT)){
    alphaT <- alpha - lambda / kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda, alphaT = alphaT,
                 r = r, time = time, ttm = ttm)

  mu.fut <- .mu.fut.schwartz2f(x0 = log(s0), delta0 = delta0, mu = mu,
                               sigmaS = sigmaS, kappa = kappa,
                               sigmaE = sigmaE, rho = rho, alpha = alpha,
                               alphaT = alphaT, r = r, time = time, ttm = ttm,
                               measure = measure)

  sigma.fut <- .sigma.fut.schwartz2f(sigmaS = sigmaS, kappa = kappa, sigmaE = sigmaE,
                                     rho = rho, time = time, ttm = ttm)

  return(dnorm(log(x), mean = mu.fut, sd = sqrt(sigma.fut), ...) / x)
}

setMethod("dfutures", signature(x = "ANY", time = "ANY", ttm = "ANY", s0 = "numeric"),
          dfutures.default)
### < ---------------------------------------------------------------------- >

dfutures.schwartz2f <- function(x, time = 0.1, ttm = 1, s0, r = 0.05, lambda = 0,
                                alphaT = NULL,
                                measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)
  
  tmp.coef <- coef(s0)

  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    alphaT <- tmp.coef$alpha
  }else if(missing(alphaT)){
    alphaT <- tmp.coef$alpha - lambda / tmp.coef$kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(dfutures.default(x, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho, r = r,
                          alphaT = alphaT, measure = measure, ...))

}

setMethod("dfutures", signature(x = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f"),
          dfutures.schwartz2f)
### < ---------------------------------------------------------------------- >

dfutures.schwartz2f.fit <- function(x, time = 0.1, ttm = 1, s0,
                                    measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)

  tmp.coef <- coef(s0)

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho
  r <- tmp.coef$r
  alphaT <- tmp.coef$alphaT

  return(dfutures.default(x, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho,
                          r = r, alphaT = alphaT, measure = measure, ...))
}

setMethod("dfutures", signature(x = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f.fit"),
          dfutures.schwartz2f.fit)
### < ---------------------------------------------------------------------- >

### < ====================================================================== >
setGeneric("qfutures",
           function(p, time = 0.1, ttm = 1, s0, ...)
           standardGeneric("qfutures"))

qfutures.default <- function(p, time = 0.1, ttm = 1, s0 = 50, delta0 = 0, mu = 0.1,
                             sigmaS = 0.3, kappa = 1, alpha = 0,
                             sigmaE = 0.5, rho = 0.75,
                             r = 0.05, lambda = 0, alphaT = NULL,
                             measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)
  
  stopifnot(time <= ttm)
  if((missing(lambda) | missing(alpha)) & missing(alphaT)){
    warning("Both 'alphaT' and ('lambda' or 'alpha') are missing!\n",
            "The mean-level of the convenience yield is set to zero.")
    alphaT <- 0
  }else if(missing(alphaT)){
    alphaT <- alpha - lambda / kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda, alphaT = alphaT,
                 r = r, time = time, ttm = ttm)

  
  mu.fut <- .mu.fut.schwartz2f(x0 = log(s0), delta0 = delta0, mu = mu,
                               sigmaS = sigmaS, kappa = kappa,
                               sigmaE = sigmaE, rho = rho, alpha = alpha,
                               alphaT = alphaT, r = r, time = time, ttm = ttm,
                               measure = measure)

  sigma.fut <- .sigma.fut.schwartz2f(sigmaS = sigmaS, kappa = kappa, sigmaE = sigmaE,
                                     rho = rho, time = time, ttm = ttm)

  return(exp(qnorm(p = p, mean = mu.fut, sd = sqrt(sigma.fut), ...)))
}

setMethod("qfutures", signature(p = "ANY", time = "ANY", ttm = "ANY", s0 = "numeric"),
          qfutures.default)
### < ---------------------------------------------------------------------- >

qfutures.schwartz2f <- function(p, time = 0.1, ttm = 1, s0, r = 0.05, lambda = 0,
                                alphaT = NULL,
                                measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)
  
  tmp.coef <- coef(s0)

  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    alphaT <- tmp.coef$alpha
  }else if(missing(alphaT)){
    alphaT <- tmp.coef$alpha - lambda / tmp.coef$kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(qfutures.default(p, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho,
                          r = r, alphaT = alphaT, measure = measure, ...))


}
setMethod("qfutures", signature(p = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f"),
          qfutures.schwartz2f)
### < ---------------------------------------------------------------------- >

qfutures.schwartz2f.fit <- function(p, time = 0.1, ttm = 1, s0,
                                    measure = c("P", "Q"), ...)
{
  measure <- match.arg(measure)

  tmp.coef <- coef(s0)

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho
  r <- tmp.coef$r
  alphaT <- tmp.coef$alphaT

  return(qfutures.default(p, time = time, ttm = ttm, s0 = tmp.coef$s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho,
                          r = r, alphaT = alphaT, measure = measure, ...))
}

setMethod("qfutures", signature(p = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f.fit"),
          qfutures.schwartz2f.fit)
### < ---------------------------------------------------------------------- >


### < ====================================================================== >
setGeneric("rfutures",
           function(n, time = 0.1, ttm = 1, s0, ...)
           standardGeneric("rfutures"))

rfutures.default <- function(n, time = 0.1, ttm = 1, s0 = 50, delta0 = 0, mu = 0.1,
                             sigmaS = 0.3, kappa = 1, alpha = 0,
                             sigmaE = 0.5, rho = 0.75, r = 0.05,
                             lambda = 0, alphaT = NULL,
                             measure = c("P", "Q"))
{
  measure <- match.arg(measure)

  stopifnot(time <= ttm)  
  if((missing(lambda) | missing(alpha)) & missing(alphaT)){
    warning("Both 'alphaT' and ('lambda' or 'alpha') are missing!\n",
            "The mean-level of the convenience yield is set to zero.")
    alphaT <- 0
  }else if(missing(alphaT)){
    alphaT <- alpha - lambda / kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda, alphaT = alphaT,
                 r = r, time = time, ttm = ttm)

  mu.fut <- .mu.fut.schwartz2f(x0 = log(s0), delta0 = delta0, mu = mu, sigmaS = sigmaS,
                               kappa = kappa, sigmaE = sigmaE, rho = rho, alpha = alpha,
                               alphaT = alphaT, r = r, time = time, ttm = ttm,
                               measure = measure)

  sigma.fut <- .sigma.fut.schwartz2f(sigmaS = sigmaS, kappa = kappa, sigmaE = sigmaE,
                                     rho = rho, time = time, ttm = ttm)

  return(exp(rnorm(n = n, mean = mu.fut, sd = sqrt(sigma.fut))))
}

setMethod("rfutures", signature(n = "ANY", time = "ANY", ttm = "ANY", s0 = "numeric"),
          rfutures.default)
### < ---------------------------------------------------------------------- >

rfutures.schwartz2f <- function(n, time = 0.1, ttm = 1, s0, r = 0.05, lambda = 0,
                                alphaT = NULL,
                                measure = c("P", "Q"))                                
{
  measure <- match.arg(measure)
  
  tmp.coef <- coef(s0)

  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    alphaT <- tmp.coef$alpha
  }else if(missing(alphaT)){
    alphaT <- tmp.coef$alpha - lambda / tmp.coef$kappa
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'lambda' is ignored.")
  }

  s0 <- tmp.coef$s0
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(rfutures.default(n, time = time, ttm = ttm, s0 = s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho, r = r,
                          alphaT = alphaT, measure = measure))

}

setMethod("rfutures", signature(n = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f"),
          rfutures.schwartz2f)
### < ---------------------------------------------------------------------- >

rfutures.schwartz2f.fit <- function(n, time = 0.1, ttm = 1, s0,
                                    measure = c("P", "Q"))
{
  measure <- match.arg(measure)
  
  tmp.coef <- coef(s0)

  s0 <- tmp.coef$s0
  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho
  r <- tmp.coef$r
  alphaT <- tmp.coef$alphaT

  return(rfutures.default(n, time = time, ttm = ttm, s0 = s0, delta0 = delta0,
                          mu = mu, sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                          sigmaE = sigmaE, rho = rho,
                          r = r, alphaT = alphaT, measure = measure))

}

setMethod("rfutures", signature(n = "ANY", time = "ANY", ttm = "ANY",
                                s0 = "schwartz2f.fit"),
          rfutures.schwartz2f.fit)
### < ---------------------------------------------------------------------- >


### < ====================================================================== >
setGeneric("filter.schwartz2f",
           function(data, ttm, s0, ...)
           standardGeneric("filter.schwartz2f"))

filter.schwartz2f.default <- function(data, ttm, s0 = 50, delta0 = 0,
                                      mu = 0.1, sigmaS = 0.3, kappa = 1, alpha = 0,
                                      sigmaE = 0.5, rho = 0.75, r = 0.05,
                                      lambda = 0, alphaT = NULL, deltat = 1/260,
                                      meas.sd = rep(1e-3, ncol(data)),
                                      P0 = 0.5 * diag(c(log(s0), abs(delta0))))
{
  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    lambda <- 0
  }else if(missing(lambda)){
    lambda <- kappa * (alpha - alphaT)
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'alphaT' is ignored.")
  }

  ## Parameters must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda,
                 r = r)

  data <- log(as.matrix(data))

  stateSpace <- .state.space.2f(y = data, ttm = ttm, deltat = deltat,
                                x0 = log(s0), delta0 = delta0,
                                kappa = kappa, mu = mu, alpha = alpha, lambda = lambda,
                                sigmaS = sigmaS, sigmaE = sigmaE, rho = rho,
                                gg = meas.sd, r = r, d = ncol(data), n = nrow(data))

  filtered.ts <- fkf(a0 = stateSpace$a0,
                     P0 = stateSpace$P0,
                     Tt = stateSpace$Tt,
                     dt = stateSpace$dt,
                     HHt = stateSpace$HHt,
                     yt = stateSpace$yt,
                     Zt = stateSpace$Zt,
                     ct = stateSpace$ct,
                     GGt = stateSpace$GGt)

  state <- cbind(S = exp(filtered.ts$att[1,]),
                 delta = filtered.ts$att[2,])

  return(list(state = state, fkf.obj = filtered.ts))
}
setMethod("filter.schwartz2f", signature(data = "ANY", ttm = "ANY",
                                     s0 = "numeric"),
          filter.schwartz2f.default)
### < ---------------------------------------------------------------------- >

filter.schwartz2f.schwartz2f <- function(data, ttm, s0, r = 0.05,
                                         lambda = 0, alphaT = NULL,
                                         deltat = 1/260,
                                         meas.sd = rep(1e-3, ncol(data)),
                                         P0 = 0.1 * diag(2))
{
  tmp.coef <- coef(s0)

  if(missing(lambda) & missing(alphaT)){
    warning("Both 'lambda' and 'alphaT' are missing!\n",
            "The market price of risk is set to zero")
    lambda <- 0
  }else if(missing(lambda)){
    lambda <- tmp.coef$kappa * (tmp.coef$alpha - alphaT)
  }else if(!missing(alphaT) & !missing(lambda)){
    warning("Both 'alphaT' and 'lambda' were passed: 'alphaT' is ignored.")
  }

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  kappa <- tmp.coef$kappa
  alpha <- tmp.coef$alpha
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho

  return(filter.schwartz2f.default(data = data, ttm = ttm,
                                   s0 = tmp.coef$s0, delta0 = delta0, mu = mu,
                                   sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                                   sigmaE = sigmaE, rho = rho, r = r,
                                   lambda = lambda, deltat = deltat, meas.sd = meas.sd))

}
setMethod("filter.schwartz2f", signature(data = "ANY", ttm = "ANY",
                                     s0 = "schwartz2f"),
          filter.schwartz2f.schwartz2f)
### < ---------------------------------------------------------------------- >

filter.schwartz2f.schwartz2f.fit <- function(data, ttm, s0)
{
  tmp.coef <- coef(s0)

  delta0 <- tmp.coef$delta0
  mu <- tmp.coef$mu
  sigmaS <- tmp.coef$sigmaS
  alpha <- tmp.coef$alpha
  kappa <- tmp.coef$kappa
  sigmaE <- tmp.coef$sigmaE
  rho <- tmp.coef$rho
  r <- tmp.coef$r
  lambda <- tmp.coef$lambda

  return(filter.schwartz2f.default(data = data, ttm = ttm,
                                   s0 = tmp.coef$s0, delta0 = delta0, mu = mu,
                                   sigmaS = sigmaS, kappa = kappa, alpha = alpha,
                                   sigmaE = sigmaE, rho = rho, r = r,
                                   lambda = lambda, deltat = s0@deltat,
                                   meas.sd = s0@meas.sd))

}

setMethod("filter.schwartz2f", signature(data = "ANY", ttm = "ANY",
                                     s0 = "schwartz2f.fit"),
          filter.schwartz2f.schwartz2f.fit)




### < ---------------------------------------------------------------------- >
futuresplot <- function(futures, type = c("forward.curve", "ttm"), ...)
                        
{
  type <- match.arg(type)
  dates <- as.Date(rownames(futures$ttm))

  if(type == "ttm"){
    col <- rainbow(ncol(futures$ttm))
    plot(dates,futures$ttm[,1], ylim = range(futures$ttm, na.rm = TRUE), type = "l", col = col[1],
         ylab = "Time to maturity [1 day]", ...)
    for(i in 2:ncol(futures$ttm)){
      lines(dates, futures$ttm[,i], col = col[i])
    }
  }else{
    plot(dates, futures$price[,1],
         xlim = c(min(dates), max(dates) + max(futures$ttm[nrow(futures$ttm),], na.rm = TRUE)),
         ylim = range(futures$price, na.rm = TRUE), type = "l", ylab = "", ...)

    col <- rainbow(10)
    col.idx <- rep(1:length(col), length = nrow(futures$price))
    tmp.data <- cbind(futures$price, futures$ttm, data.frame(col.idx),
                      1:nrow(futures$price))

    plot.forward.curve <- function(x, col, d, dates){
      lines(dates[x[2 * d + 2]] + c(0,x[(d + 1):(2 * d)]), x[c(1,1:d)], 
            col = col[x[2 * d + 1]], lty = "dotted")
    }
    apply(tmp.data, 1, plot.forward.curve, d = ncol(futures$ttm), col = col,
          dates = dates)
    legend("topleft", legend = c("Closest to maturity contract", "Forward Curves"),
           lty = c("solid", "dotted"), bg = "white")

  }
}
