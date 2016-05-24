DiscreteLaplace <-
function() {
  make.DiscreteLaplace.link <- function() {

    # A name to be used for the link
    link <- "DiscreteLaplace"
    
    # Link function
    linkfun <- function(mu) log( (sqrt(1 + mu^2) - 1) / mu)
    
    # Inverse link function
    linkinv <- function(theta) {
      pmax(2*exp(theta) / (1 - exp(2*theta)), .Machine$double.eps)
    }
    
    # Derivative function(theta) dmu/dtheta
    mu.theta <- function(theta) {
      pmax( (2*exp(theta)*(1+exp(2*theta))) / ((exp(2*theta)-1)^2), .Machine$double.eps)
    }
    
    # TRUE if theta is in the domain of linkinv
    validtheta <- function(theta) all(theta < 0 & !is.infinite(theta))

    structure(list(linkfun = linkfun, linkinv = linkinv, mu.theta = mu.theta, 
        valideta = validtheta, name = link), class = "link-glm")
  }

  stats <- make.DiscreteLaplace.link()
  linktemp <- stats$name

  variance <- function(mu) mu * sqrt(1 + mu^2)
  validmu <- function(mu) return(mu > 0 & !is.infinite(mu))

  loglikeh <- function(mu, d) {
    p <- (sqrt(1 + mu^2) - 1) / mu
    return(log(1-p) - log(1+p) + d*log(p))    
  }

  dev.res <- function(d, mu, wt) {
    mu[mu < 10e-8] <- 10e-8
    
    muinv <- 1/mu
    
    ret.val <- rep(0, length(d))
    d.zero <- d == 0
    d.not.zero <- !d.zero
    
    if (sum(d.zero) > 0) {
      const <- muinv[d.zero] * (sqrt((mu[d.zero])^2+1) - 1)
      ret.val[d.zero] <- 2*log( (1 + const) / (1 - const) )
    }

    if (sum(d.not.zero) > 0) {
      d.loglikeh <- d[d.not.zero]
      ret.val[d.not.zero] <- 2 * (loglikeh(d.loglikeh, d.loglikeh) - loglikeh(mu[d.not.zero], d.loglikeh))
    }
    
    return(list(ret.val = ret.val, mu = mu))
  }
  
  dev.null <- function(d, wt) {    
    ret.val <- rep(0, length(d))
    d.zero <- d == 0
    d.not.zero <- !d.zero
    
    mu <- rep(mean(d), length(d), na.rm = TRUE)
    muinv <- 1/mu
    
    if (sum(d.zero) > 0) {
      const <- muinv[d.zero] * (sqrt((mu[d.zero])^2+1) - 1)
      ret.val[d.zero] <- 2*log( (1 + const) / (1 - const) )
    }

    if (sum(d.not.zero) > 0) {
      d.loglikeh <- d[d.not.zero]
      #ret.val[d.not.zero] <- - 2 * loglikeh(mu[d.not.zero], d.loglikeh)
      #ret.val[d.not.zero] <- 2 * (loglikeh(d.loglikeh, d.loglikeh) - loglikeh(mu, d.loglikeh))
      ret.val[d.not.zero] <- 2 * (loglikeh(d.loglikeh, d.loglikeh) - loglikeh(mu[d.not.zero], d.loglikeh))
    }
    
    return(list(ret.val = ret.val, mu = mu))
  }
  
  dev.resids <- function(y, mu, wt) {
    d <- y
    mu.org <- mu
  
    if (length(mu) == 1) {
      mu <- rep(mu, length(y))
    }
  
    if (length(mu) != length(d)) {
      print(mu)
      stop()
    }

    if (any(d < 0)) stop("any(d < 0)")
    
    ret <- NULL
    
    if (all(is.infinite(mu))) {
      ret <- dev.null(d, wt)
    } else {
      ret <- dev.res(d, mu, wt)
    }
    
    ret.val <- ret$ret.val
    ret.val.wt <- wt * ret.val

    return(ret.val.wt)
  }
  
  aic <- function(y, n, mu, wt, dev) {
    return(-2 * sum(loglikeh(mu, y) * wt))
  }

  initialise <- expression({
    if (any(y < 0)) stop("negative values of d not allowed")    
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })

  structure(list(family = "DiscreteLaplace", link = linktemp, linkfun = stats$linkfun, 
      linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, aic = aic,
      mu.eta = stats$mu.theta, initialize = initialise, validmu = validmu, valideta = stats$valideta), 
      class = "family")
}

