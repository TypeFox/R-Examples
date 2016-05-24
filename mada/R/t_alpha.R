logit<-function(p) 
{
    out <- p
    inds <- seq_along(p)[!is.na(p)]
    if (any((p[inds] < 0) | (p[inds] > 1))) 
        stop("invalid proportions input")
    out[inds] <- log(p[inds]/(1 - p[inds]))
    out[inds][p[inds] == 0] <- -Inf
    out[inds][p[inds] == 1] <- Inf
    out
}

inv.logit<-function(x) plogis(x)

### t_alpha as link-function

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
trafo <- function(alpha, x){return(talpha(alpha)$linkfun(x))}
inv.trafo <- function(alpha, x){return(talpha(alpha)$linkinv(x))}
jacobiantrafo <- function(alpha,x)return(alpha/x + (2-alpha)/(1-x))
