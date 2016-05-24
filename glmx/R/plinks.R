## angular transformation (aka arcsine transformation)

angular <- function(verbose = FALSE)
{
  linkfun <- function(mu) asin(sqrt(mu))  
  linkinv <- function(eta) {
    etastar <- pmin(asin(1) - .Machine$double.eps, pmax(.Machine$double.eps, eta))
    if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
      warning("truncation in inverse link function")
    return((sin(etastar))^2)
  }      
  mu.eta <- function(eta) 2 * cos(eta) * sin(eta)  
  valideta <- function(eta) {
    if(verbose && (!all(eta <= asin(1) & 0 <= eta)))
      warning("some of the current etas are out of range")
    TRUE
  }
  name <- "angular"
 
  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}



## just for completeness: log-log link in addition to complementary log-log

loglog <- function() {
  structure(list(
    linkfun = function(mu) -log(-log(mu)),
    linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
    mu.eta = function(eta) {
      eta <- pmin(eta, 700)
      pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
    },
    valideta = function(eta) TRUE,
    name = "loglog"),
    class = "link-glm")
}

## transformations from Aranda-Ordaz (1981)
## ao1 is the symmetric link from Section 2.1

ao1 <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(verbose && phi < 0) warning("sign of phi ignored as ao1(phi) = ao1(-phi)")
  phi <- abs(phi)
  if(phi == 0) {
    rval <- make.link("logit")
    rval$name <- "ao1"
    return(rval)
  }

  linkfun <- function(mu) (2/phi) * (mu^phi - (1 - mu)^phi)/(mu^phi + (1 - mu)^phi)
  
  linkinv <- function(eta) {
    etastar <- pmin(2/phi - .Machine$double.eps, pmax(-2/phi + .Machine$double.eps, eta))
    if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
      warning("truncation in inverse link function")
    ((1 + 0.5 * phi * etastar)^(1/phi))/((1 + 0.5 * phi * etastar)^(1/phi) + (1 - 0.5 * phi * etastar)^(1/phi))
  }
      
  mu.eta <- function(eta) {
    phieta1 <- 1 + 0.5 * phi * eta
    phieta2 <- 1 - 0.5 * phi * eta
    phieta1^((1/phi) - 1) * 0.5/(phieta1^(1/phi) + phieta2^(1/phi)) - (phieta1^(1/phi)) *
      (phieta1^((1/phi) - 1) * 0.5 - phieta2^((1/phi) - 1) * (0.5))/(phieta1^(1/phi) + phieta2^(1/phi))^2
  }
  
  valideta <- function(eta) {
    if(verbose && !all(abs(0.5 * phi * eta) < 1)) warning("some of the current etas are out of range")
    TRUE
  }
  
  name <- "ao1" ## "Aranda-Ordaz symmetric"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}
    
## transformations from Aranda-Ordaz (1981)
## ao2 is the assymmetric link from Section 2.2

ao2 <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(phi == 1) {
    rval <- make.link("logit")
    rval$name <- "ao2"
    return(rval)
  }
  if(phi == 0) {
    rval <- make.link("cloglog")
    rval$name <- "ao2"
    return(rval)
  }

  linkfun <- function(mu) log(((1 - mu)^(-phi) - 1)/phi)
    
  linkinv <- function(eta){
    if(phi < 0) {
      etastar <- pmin(eta, log(-1/phi) - .Machine$double.eps)
      if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
        warning("truncation in inverse link function")
      eta <- etastar
    }
    1- (1 + phi * exp(eta))^(-1/phi)
  }
      
  mu.eta <- function(eta) exp(eta) * (1 + phi * exp(eta))^(-(1 + phi)/phi)

  valideta <- function(eta) {
    if(verbose && !all(phi * exp(eta) > -1)) warning("some of the current etas are out of range")
    TRUE
  }
  
  name <- "ao2" ## "Aranda-Ordaz asymmetric"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}



## folded exponential transformation (Piepho 2003)

foldexp <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(phi == 0) {
    structure(list(
      linkfun = function(mu) mu - 0.5,
      linkinv = function(eta) {
        etastar <- pmin(0.5 - .Machine$double.eps, pmax(-0.5 + .Machine$double.eps, eta))
        if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
          warning("truncation in inverse link function")
        etastar + 0.5
      },
      mu.eta = function(eta) rep(1, length(eta)),
      valideta = function(eta) TRUE,
      name = "foldexp"),
      class = "link-glm")  
  }
  
  upper <- (exp(phi) - 1)/(2 * phi)
  lower <- (1 - exp(phi))/(2 * phi)

  linkfun <- function(mu) (exp(phi * mu) - exp(phi * (1 - mu)))/(2 * phi)
      
  dfoldexp <- function(mu) (exp(phi * mu) + exp(phi * (1 - mu)))/2
  
  linkinv <- function(eta) {
    etastar <- pmin(upper - .Machine$double.eps, pmax(lower + .Machine$double.eps, eta))
    if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
      warning("truncation in inverse link function")
    0.5 + asinh((phi * etastar)/exp(phi/2))/phi
  }
      
  mu.eta <- function(eta) 1/dfoldexp(linkinv(eta))
  
  valideta <- function(eta) TRUE

  name <- "foldexp" ## "folded exponential"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}

## Guerrero-Johnson (1982) transformation 

gj <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(phi == 0) {
    rval <- make.link("logit")
    rval$name <- "gj"
    return(rval)
  }

  linkfun <- function(mu) ((mu/(1 - mu))^phi - 1)/phi
  
  linkinv <- function(eta){
    etastar <- pmax((-sign(phi)/phi) + .Machine$double.eps, eta)
    if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
      warning("truncation in inverse link function")
    ((phi * etastar + 1)^(1/phi))/(1 + ((phi * etastar + 1)^(1/phi)))
  }
      
  mu.eta <- function(eta) {
    pe1 <- (phi * eta + 1)^((1/phi) - 1)
    pe2 <- (phi * eta + 1)^(1/phi)
    pe1/(1 + pe2) - pe2 * (pe1/(1 + pe2)^2)
  }
  
  valideta <- function(eta) TRUE
    ## old version:
    ##  all(-sign(phi)/phi < eta)

  name <- "gj" ## "Guerrero-Johnson"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}


## t_alpha link function

talpha <- function(alpha, verbose = FALSE,
  splineinv = TRUE, eps = 2 * .Machine$double.eps, maxit = 100)
{
  ## parameter processing
  if(alpha > 2 | alpha < 0) stop("alpha must be between 0 and 2")
  if(alpha == 1) {
    rval <- make.link("logit")
    rval$name <- "talpha"
    return(rval)  
  }

  linkfun <- function(mu) alpha * log(mu) - (2 - alpha) * log(1 - mu)
        
  dtrafo <- function(alpha, x) alpha/x + (2 - alpha)/(1 - x)
  
  linkinv <- function(eta) {
    ## inverse link is easy to compute for alpha 0, 1, 2
    ## alpha == 1 already special-cased above
    if(alpha == 0){
      etastar <- pmax(eta, .Machine$double.eps)
      if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
        warning("truncation in inverse link function")
      return(1 - exp(-etastar/2))
    }
    if(alpha == 2) {
      etastar <- pmin(eta, -.Machine$double.eps)
      if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
        warning("truncation in inverse link function")
      return(exp(etastar/2))
    }

    ## otherwise some numeric technique is required:
    ## either via splines (quick & dirty) or via Newton (cleaner but somewhat slower)
    if(splineinv) {

      ## this is somewhat less precise than the Newton solution
      ## but faster for moderately to large data sets
      mu0 <- plogis(seq(-30, 30, by = 0.01))
      eta0 <- linkfun(mu0)
      mu <- splinefun(eta0, mu0)(eta)

    } else {

      ## simple newton algorithm (maxit is rarely reached, only when alpha close to 0 or 2)
      ## 8 times faster than uniroot
      mu <- plogis(eta)
      ok <- FALSE
      n <- 0
      while(!ok) {
        mu.prev <- mu
        mu <- pmax(pmin(mu - (linkfun(mu) - eta)/(dtrafo(alpha, mu) - 1), 1 - eps), eps)
        ok <- max(abs(mu.prev - mu)) < eps
        if(n > maxit) {
          warning("maxit reached")
	  ok <- TRUE
        }
        n <- n + 1
      }

    }
    mu
  }
  
  mu.eta <- function(eta) 1/dtrafo(alpha, linkinv(eta))

  valideta <- function(eta) TRUE
  ## old version:
  ##    if(alpha == 0){  if(any(eta <= 0)){FALSE}
  ##            if(any(eta > 0)){TRUE}}
  ##    if(alpha == 2){  if(any(eta >= 0)){FALSE}
  ##            if(any(eta < 0)){TRUE}}
  
  name <- "talpha"
  
  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}


## Rocke (1993) beta transformation as link function

rocke <- function(shape1, shape2, verbose = FALSE)
{
  ## parameter processing
  if(missing(shape2)) {
    shape2 <- rep(shape1, length.out = 2)[2]
    shape1 <- shape1[1]
  }  
  if(shape1 == 0 & shape2 == 0) {
    rval <- make.link("logit")
    rval$name <- "rocke"
    return(rval)
  }
  if(shape1 < 0 | shape2 < 0) stop("parameters of Rocke transformation need to be both strictly positive or both equal to 0")
  # what to do if shape1 == 0, and shape2 > 0? nothing: then the integral does not exist!
    
  shift <- pbeta(0.5, shape1, shape2)
  betafactor <- beta(shape1, shape2)
      
  linkfun <- function(mu) betafactor * (pbeta(mu, shape1, shape2) - shift)
  
  linkinv <- function(eta){
    etastar <- pmin((1-shift)*betafactor, pmax(eta, -shift*betafactor))  
    if(verbose && !isTRUE(all.equal(as.vector(eta), as.vector(etastar))))
      warning("truncation in inverse link function")
    pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps,
      qbeta(etastar/betafactor + shift, shape1, shape2)))
  }
  
  mu.eta <- function(eta)
    1/(betafactor * (dbeta(linkinv(eta), shape1, shape2)))

  valideta <- function(eta) {
    if(verbose && any( eta < -shift * betafactor | eta > (1 - shift) * betafactor ))
      warning("some of the current etas are out of range")
    TRUE
  }
  
  name <- "rocke" ## "Rocke's beta"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}


## some form of generalized logistic link
## constructed as right-censored negative
## binomial distribution with parameter theta

nblogit <- function(theta) {
  ## parameter processing
  if(theta <= 0) stop("theta must be positive")
  if(theta == 1) {
    rval <- make.link("logit")
    rval$name <- "nblogit"
    return(rval)
  }
  if(!is.finite(theta)) {
    rval <- make.link("cloglog")
    rval$name <- "nblogit"
    return(rval)
  }

  linkfun <- function(mu) log(1 - (1 - mu)^(1/theta)) - log(1 - mu)/theta + log(theta)
  linkinv <- function(eta) 1 - (1 / (1 + exp(eta)/theta))^theta
  mu.eta <- function(eta) exp( -(theta + 1) * log(1 + exp(eta)/theta) + eta )
  valideta <- function(eta) TRUE
  name <- "nblogit"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}

gosset <- function(nu, verbose = FALSE)
{
  if(nu <= 0) stop("nu must be positive")
  if(verbose && nu < 0.2) warning("gosset link might be unreliable for nu < 0.2")
  qqt <- function(p, nu) sign(p - 0.5) * sqrt(qf(1 - 2 * pmin(p, 1 - p), 1, nu))

  linkfun <- function(mu) qqt(mu, nu)
  linkinv <- function(eta) {
    thresh <- -qqt(.Machine$double.eps, nu)
    eta <- pmin(thresh, pmax(eta, -thresh))
    pt(eta, nu)
  }
  mu.eta <- function(eta) pmax(dt(eta, nu), .Machine$double.eps)
  valideta <- function(eta) TRUE
  name <- "gosset"
  
  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}
       
pregibon <- function(a, b) {

  ## allow specification via vector as first argument
  if(missing(b)) {
    b <- rep(a, length.out = 2)[2]
    a <- a[1]
  }

  linkfun <- function(mu) -qpregibon(1 - mu, a = a, b = b)

  linkinv <- function(eta) {
    eps <- .Machine$double.eps^0.5
    tlo <- qpregibon(eps, a = a, b = b)
    thi <- qpregibon(1 - eps, a = a, b = b)
    eta <- -pmin(thi, pmax(-eta, tlo))
    1 - ppregibon(-eta, a = a, b = b)
  }

  mu.eta <- function(eta) pmax(dpregibon(-eta, a = a, b = b), .Machine$double.eps^0.5)
  valideta <- function(eta) TRUE
  name <- "pregibon"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}
