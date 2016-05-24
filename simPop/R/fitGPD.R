# Taken from the now archived R package POT
#Descripton of the last available version:
#Package: POT
#Version: 1.1-3
#Date: 2012-10-30
#Title: Generalized Pareto Distribution and Peaks Over Threshold
#Author: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
#    Maintainer: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
#    Depends: R (>= 1.8.0)
#Description: Some functions useful to perform a Peak Over Threshold
#analysis in univariate and bivariate cases. A user's guide is
#    available.
#    License: GPL (>= 2)
#    URL: http://r-forge.r-project.org/projects/pot/
#    Repository: CRAN
#    Repository/R-Forge/Project: pot
#    Repository/R-Forge/Revision: 492
#    Repository/R-Forge/DateTimeStamp: 2012-10-30 14:21:03
#    Date/Publication: 2012-11-06 09:49:26
#    Packaged: 2012-10-30 15:22:42 UTC; rforge
    

## In this file, several functions to estimates the GPD parameters
## are available:
##   1) Moments Estimator
##   2) Unbiased Probability Weighted Moment (PWMU) Estimator
##   3) Biased Probability Weighted Moment (PWMB) Estimator
##   4) Maximum Likelihood Estimator
##   5) Pickands' Estimator
##   6) Minimum Density Power Divergence Estimator
##   7) Method of Medians Estimator
##   8) Likelihood Moment Estimator
##   9) Maximum Goodness-of-Fit Estimator

## A generic function for estimate the GPD parameters
fitgpd <- function(data, threshold, est = "mle", ...){
  threshold.call <- deparse(threshold)
  if (!(est %in% c("moments", "pwmb", "pwmu", "mle", "pickands",
                   "mdpd", "med", "lme", "mgf", "mple")))
    stop("Unknown estimator. Please check the ``est'' argument.")
  
  fitted <- switch(est, 'moments' = gpdmoments(data, threshold),
                   'pwmb' = gpdpwmb(data, threshold, ...),
                   'pwmu' = gpdpwmu(data, threshold),
                   'mle' = gpdmle(data, threshold, ...),
                   'pickands' = gpdpickands(data, threshold),
                   'mdpd' = gpdmdpd(data, threshold, ...),
                   'med' = gpdmed(data, threshold, ...),
                   'lme' = gpdlme(data, threshold, ...),
                   'mgf' = gpdmgf(data, threshold, ...),
                   'mple' = gpdmple(data, threshold, ...)
                   )
  fitted$threshold.call <- threshold.call
  class(fitted) <- c("uvpot","pot")
  return(fitted)
}

##Maximum penalized likelihood estimator
gpdmple <- function(x, threshold, start, ..., std.err.type =
                    "observed", corr = FALSE, method = "BFGS",
                    warn.inf = TRUE, alpha = 1, lambda = 1){

  if (all(c("observed", "expected", "none") != std.err.type))
    stop("``std.err.type'' must be one of 'observed', 'expected' or 'none'")
  
  nlpot <- function(scale, shape) { 
    ans <- -.C("gpdlik", exceed, nat, threshold, scale,
                shape, dns = double(1), PACKAGE = "simPop")$dns

    if ((shape > 0) & (shape <1))
     ans <- lambda * ((1 / (1 - shape) - 1)^alpha) + ans

    if (shape >= 1)
      ans <- 1e6

    return(ans)
  }
  
  nn <- length(x)
  
  threshold <- rep(threshold, length.out = nn)
  
  high <- (x > threshold) & !is.na(x)
  threshold <- as.double(threshold[high])
  exceed <- as.double(x[high])
  nat <- length(exceed)
  
  if(!nat) stop("no data above threshold")
  
  pat <- nat/nn
  param <- c("scale", "shape")
  
  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0)
    start$scale <- mean(exceed) - min(threshold)
    
    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nlpot)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlpot) <- c(f[m], f[-m])
  nllh <- function(p, ...) nlpot(p, ...)
  
  if(l > 1)
    body(nllh) <- parse(text = paste("nlpot(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if( warn.inf && do.call("nllh", start.arg) == 1e6 )
    warning("negative log-likelihood is infinite at starting values")
  
  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
  
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  if (std.err.type != "none"){
    
    tol <- .Machine$double.eps^0.5
    
    if(std.err.type == "observed") {
      
      var.cov <- qr(opt$hessian, tol = tol)
      if(var.cov$rank != ncol(var.cov$qr)){
        warning("observed information matrix is singular; passing std.err.type to ``expected''")
        obs.fish <- FALSE
        return
      }
      
      if (std.err.type == "observed"){
        var.cov <- solve(var.cov, tol = tol)
        
        std.err <- diag(var.cov)
        if(any(std.err <= 0)){
          warning("observed information matrix is singular; passing std.err.type to ``expected''")
          std.err.type <- "expected"
          return
        }
        
        std.err <- sqrt(std.err)
        
        if(corr) {
          .mat <- diag(1/std.err, nrow = length(std.err))
          corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
          diag(corr.mat) <- rep(1, length(std.err))
        }
        else {
          corr.mat <- NULL
        }
      }
    }
    
    if (std.err.type == "expected"){
      
      shape <- opt$par[2]
      scale <- opt$par[1]
      a22 <- 2/((1+shape)*(1+2*shape))
      a12 <- 1/(scale*(1+shape)*(1+2*shape))
      a11 <- 1/((scale^2)*(1+2*shape))
      ##Expected Matix of Information of Fisher
      expFisher <- nat * matrix(c(a11,a12,a12,a22),nrow=2)

      expFisher <- qr(expFisher, tol = tol)
      var.cov <- solve(expFisher, tol = tol)
      std.err <- sqrt(diag(var.cov))
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      else
        corr.mat <- NULL
    }

    colnames(var.cov) <- nm
    rownames(var.cov) <- nm
    names(std.err) <- nm
  }

  else{
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }
  
  
  param <- c(opt$par, unlist(fixed.param))
  scale <- param["scale"]
  
  var.thresh <- !all(threshold == threshold[1])

  if (!var.thresh)
    threshold <- threshold[1]
  
  list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = unlist(fixed.param), param = param,
       deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
       counts = opt$counts, message = opt$message, threshold = threshold,
       nat = nat, pat = pat, data = x, exceed = exceed, scale = scale,
       var.thresh = var.thresh, est = "MPLE", logLik = -opt$value,
       opt.value = opt$value)
}


## Maximum goodness-of-fit estimator
gpdmgf <- function(x, threshold, start, stat, ...,
                   method = "BFGS", warn.inf = TRUE){

  nn <- length(x)
  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])

  nat <- length(exceed)
  if (!nat) 
    stop("no data above threshold")

  if (!(stat %in% c("KS","CM","AD","ADR","ADL","AD2R","AD2L",
                   "AD2")))
    stop("`stat' must be one of 'KS','CM','AD','ADR','ADL','AD2R','AD2L', 'AD2'.")

  pat <- nat/nn
  param <- c("scale", "shape")

  excess <- exceed - threshold

  excess <- sort(excess)

  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0)
    start$scale <- mean(exceed) - min(threshold)
    
    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  if (stat == "KS")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else
        1 / 2 / nat + max(abs(pgpd(excess, 0, scale, shape) -
                              ppoints(nat)))
    }
    
  if (stat == "CM")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else
        1/12/nat + sum((pgpd(excess, 0, scale, shape) -
                        ppoints(nat))^2)
    }

  if (stat == "AD")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if ((pgpd(max(excess),  0, scale, shape) >= 1) ||
            (pgpd(min(excess),  0, scale, shape) <= 0))
         1e6

        else
          -nat - mean((2*1:nat - 1) *
                      (log(pgpd(excess, 0, scale, shape)) +
                       log(1 - pgpd(rev(excess), 0, scale, shape))))
      }
    }

  if (stat == "ADR")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(max(excess),  0, scale, shape) >= 1)
         1e6

        else
          nat / 2 - 2 * sum(pgpd(excess, 0, scale, shape)) -
            mean((2 * 1:nat -1)*log(1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }
  
  if (stat == "ADL")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(min(excess),  0, scale, shape) <= 0)
         1e6

        else
          - 3 * nat / 2 + 2 * sum(pgpd(excess, 0, scale, shape)) -
            mean((2 * 1:nat -1)*log(pgpd(excess, 0, scale, shape)))
      }
    }

  if (stat == "AD2R")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(max(excess),  0, scale, shape) >= 1)
         1e6

        else
          2 * sum(log(1 - pgpd(excess, 0, scale, shape))) +
            mean((2* 1:nat - 1) / (1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }

  if (stat == "AD2L")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(min(excess),  0, scale, shape) <= 0)
         1e6

        else
          2 * sum(log(pgpd(excess, 0, scale, shape))) +
            mean((2 * 1:nat - 1) / pgpd(excess, 0, scale, shape))
      }
    }

  if (stat == "AD2")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if ((pgpd(max(excess),  0, scale, shape) >= 1) ||
            (pgpd(min(excess),  0, scale, shape) <= 0))
         1e6

        else
          2 * sum(log(pgpd(excess, 0, scale, shape)) +
                  log(1 - pgpd(excess, 0, scale, shape))) +
                    mean((2 * 1:nat - 1) / pgpd(excess, 0, scale, shape) +
                         (2 * 1:nat - 1) /
                         (1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }

  nm <- names(start)
  l <- length(nm)
  f <- formals(fun)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(fun) <- c(f[m], f[-m])
  mgf <- function(p, ...) fun(p, ...)
  
  if(l > 1)
    body(mgf) <- parse(text = paste("fun(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if( warn.inf && do.call("mgf", start.arg) == 1e6 )
    warning("Maximum goodness-of-fit function is infinite at starting values")
  
  opt <- optim(start, mgf, hessian = TRUE, ..., method = method)
  
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  tol <- .Machine$double.eps^0.5
    
  param <- c(opt$par, unlist(fixed.param))
  scale <- param["scale"]
  
  var.thresh <- !all(threshold == threshold[1])

  if (!var.thresh)
    threshold <- threshold[1]

  std.err <- std.err.type <- corr.mat <- NULL
  var.cov <- NULL

  list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = unlist(fixed.param), param = param,
       corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
       message = opt$message, threshold = threshold, nat = nat, pat = pat,
       data = x, exceed = exceed, scale = scale, var.thresh = var.thresh,
       est = "MGF", opt.value = opt$value, stat = stat)
}

##Likelihood moment estimation
gpdlme <- function(x, threshold, r = -.5, start, ...,
                   method = "BFGS"){

  nn <- length(x)
  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])

  nat <- length(exceed)
  if (!nat) 
    stop("no data above threshold")

  pat <- nat/nn
  
  excess <- exceed - threshold
  fun <- function(x){
    if (x >= 1/max(excess))
      return(1e6)
    p <- r / mean(log(1 - x * excess))
    abs(mean((1 - x * excess)^p) - 1 / (1 - r))
  }

   if (missing(start))
    start <- list(x = -1)
   
  opt <- optim(start, fun, hessian = TRUE, ..., method = method)

  if (opt$convergence != 0){
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  counts <- opt$counts
  b <- opt$par
  zero <- opt$value
  
  shape <- mean(log(1 - b*excess))
  scale <- - shape / b

  param <- c(scale, shape)
  names(param) <- c("scale", "shape")

  a11 <- scale^2 * (2 + ((r - shape)^2 + 2 * shape) /
                    (1 - 2 * r))
  a12 <- scale * (1 + (r^2 + shape^2 + shape) /
                  (1 - 2 * r))
  a22 <- (1 - r) * (1 + (2*shape^2 - 2 * shape + r) /
                    (1 - 2 * r))

  var.cov <- matrix(c(a11, a12, a12, a22), 2) / nat
  colnames(var.cov) <- c("scale", "shape")
  rownames(var.cov) <- c("scale", "shape")

  std.err <- sqrt(diag(var.cov))

  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c("scale", "shape")
  rownames(corr) <- c("scale", "shape")

  if (shape < -0.5) 
        message <- "Assymptotic theory assumptions\nfor standard error may not be fullfilled !"
    else message <- NULL
  
  var.thresh <- FALSE
  return(list(fitted.values = param, std.err = std.err, std.err.type = "expected",
              var.cov = var.cov, param = param, message = message, data = x,
              threshold = threshold, corr = corr, convergence = opt$convergence,
              counts = counts, nat = nat, pat = pat, exceed = exceed, scale = scale,
              var.thresh = var.thresh, est = "LME", opt.value = opt$value))
}

##Pickand's Estimator
gpdpickands <- function(data, threshold){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'pickands'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  excess <- sort(exceed - threshold)
  
  n <- length(excess)
  xn.2 <- excess[ceiling(n/2)]
  x3n.4 <- excess[ceiling(.75*n)]
  d <- xn.2^2 / (2 * xn.2 - x3n.4)
  
  shape <- -log(xn.2 / (x3n.4 - xn.2) ) / log(2)
  scale <- -shape * d
  
  if ( (max(excess) * shape) > -scale)
    message <- "Estimates are valid"
  
  else
    message <- "Estimates are not valid"
  
  estim <- param <- c(scale = scale, shape = shape)
  std.err <- var.cov <- corr <- NULL
  convergence <- counts <- NA
  var.thresh <- FALSE
  
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              nat = nat, pat = pat, convergence = convergence,
              corr = corr, counts = counts, exceed = exceed,
              scale = scale, var.thresh = var.thresh, est = "pickands"))
}

## Moments Estimator

gpdmoments <- function(data, threshold){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'moments'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( nat == 0 )
    stop("None observation above the specified threshold !!!")
  
  exceed <- sort(exceed)
  
  loc <- threshold
  
  ## Evaluate the excess above the threshold 
  exces <- exceed - loc
  
  m <- mean(exces)
  v <- var(exces)
  
  scale <- m / 2 * ( m^2 / v +1 )
  shape <- - ( m^2 / v -1 ) / 2
  
  estim <- param <- c(scale  = scale, shape = shape)
  convergence <- counts <- NA
  
  a11 <- 2*scale^2 * ( 1 - 6*shape + 12*shape^2)
  a12 <- - scale * (1-2*shape) * (1-4*shape+12*shape^2)
  a21 <- a12
  a22 <- (1-2*shape)^2 * (1-shape+6*shape^2)
  
  var.cov <- (1 - shape)^2 / ( (1-2*shape)*(1-3*shape)*(1-4*shape)*nat ) *
    matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)                    
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.25 ) message <- 'Assymptotic theory assumptions
for standard error may not be fullfilled !'
  else message <- NULL
  
  var.thresh <- FALSE
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              nat = nat, pat = pat, convergence = convergence,
              corr= corr, counts = counts, exceed = exceed,
              scale=scale, var.thresh = var.thresh, est = "moments"))
}

##PWMB Estimator

gpdpwmb <- function(data, threshold, a=0.35, b=0, hybrid = FALSE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'pwmb'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( nat == 0 )
    stop("None observation above the specified threshold !!!")
  
  exceed <- sort(exceed)
  
  loc <- threshold
  
  excess <- exceed - loc
  
  m <- mean(excess)
  n <- length(excess)
  p <- (1:n - a) / (n + b)
  
  t <- sum((1-p)*excess)/n
  
  shape <- - m / (m- 2*t ) + 2
  scale <- 2 * m * t / (m - 2*t )
  est <- 'PWMB'

  if (hybrid)
    if ( (max(excess) >= (-scale / shape)) & (shape < 0) ){
      shape <- -scale / max(excess)
      est <- 'PWMB Hybrid'
    }
  
  estim <- c(scale  = scale, shape = shape)
  param <-  c(scale = scale, shape = shape)
  convergence <- NA
  counts <- NA
  
  a11 <- scale^2 * (7-18*shape+11*shape^2-2*shape^3)
  a12 <- - scale * (2-shape) * (2-6*shape+7*shape^2-2*shape^3)
  a21 <- a12
  a22 <- (1-shape) * (2 -shape)^2 * (1-shape+2*shape^2)
  
  var.cov <- 1 / ( (1-2*shape) * (3-2*shape)*nat ) *
    matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.5 )
    message <- "Assymptotic theory assumptions for standard error may not be fullfilled !"
  else message <- NULL
  
  var.thresh <- FALSE
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              corr = corr, convergence = convergence, counts = counts,
              nat = nat, pat = pat, exceed = exceed,
              scale=scale, var.thresh = var.thresh, est = est))
}


## PWMU Estimator
## First, we need a function which computes the samples L-moments

samlmu <- function (x, nmom = 4, sort.data = TRUE)
{
  xok <- x[!is.na(x)]
  n <- length(xok)
  if (nmom <= 0) return(numeric(0))
  if (nmom <= 2) rnames <- paste("l", 1:nmom, sep = "_")
  else rnames <- c("l_1", "l_2", paste("t", 3:nmom, sep = "_"))
  lmom <- rep(NA, nmom)
  names(lmom) <- rnames
  if (n == 0) return(lmom)
  if (sort.data == TRUE) xok <- sort(xok)
  nmom.actual <- min(nmom, n)
  
  lmom <- .C("samlmu", as.double(xok), as.integer(nmom.actual),
             as.integer(n), lmom = double(nmom.actual),
             PACKAGE = "simPop")$lmom
  names(lmom) <- rnames
  return(lmom)
}

gpdpwmu <- function(data,threshold, hybrid = FALSE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'pwmu'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  
  if ( length(exceed) == 0 )
    stop("None observation above the specified threshold !!!")
  
  exceed <- sort(exceed)
  nat <- length( exceed )
  pat <- nat / length( data )
  
  loc <- threshold
  excess <- exceed - loc
  
  lmoments <- samlmu(excess, nmom=2, sort.data = FALSE)
  shape <- - lmoments[1]/lmoments[2] + 2
  scale <- (1 - shape)*lmoments[1]
  names(shape) <- NULL
  names(scale) <- NULL
  est <- "PWMU"
  
  if (hybrid)
    if ( (excess[nat] >= (-scale / shape)) & (shape < 0) ){
      shape <- -scale / excess[nat]
      est <- 'PWMU Hybrid'
    }
  
  
  estim <- param <- c(scale  = scale, shape = shape)
  convergence <- counts <- NA
  
  a11 <- scale^2 * (7-18*shape+11*shape^2-2*shape^3)
  a12 <- - scale * (2-shape) * (2-6*shape+7*shape^2-2*shape^3)
  a21 <- a12
  a22 <- (1-shape) * (2 -shape)^2 * (1-shape+2*shape^2)
  
  var.cov <- 1 / ( (1-2*shape) * (3-2*shape)*nat ) * matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)                    
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.5 ) message <- "Assymptotic theory assumptions
for standard error may not be fullfilled !"
  else message <- NULL
  
  var.thresh <- FALSE
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              corr = corr, convergence = convergence, counts = counts,
              nat = nat, pat = pat, exceed = exceed,
              scale=scale, var.thresh = var.thresh, est = est))
}

##MDPD estimators for the GPD.
gpdmdpd <- function(x, threshold, a, start, ...,
                    method = "BFGS", warn.inf = TRUE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'mdpd'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  if (missing(a))
    a <- .1
  
  nn <- length(x)
  
  threshold <- rep(threshold, length.out = nn)
  
  high <- (x > threshold) & !is.na(x)
  threshold <- as.double(threshold[high])
  exceed <- as.double(x[high])
  nat <- length(exceed)
  
  excess <- exceed - threshold
  
  if(!nat) stop("no data above threshold")
  
  pat <- nat/nn
    
  if(missing(start)) {
    start <- list(scale = 0, shape = 0.01)
    start$scale <- mean(exceed) - min(threshold)
  }

  start <- c(scale = start$scale, shape = start$shape)
  
  pddf <- function(param){
    ## Evaluates the (P)ower (D)ensity (D)ivergence (F)unction which is
    ## criterion function of the MDPDE
    scale <- param[1]
    shape <- param[2]
    
    if ( ((max(excess)  < (-scale / shape)) && (shape < 0)) ||
        (shape > 0) ){
      y <- pmax(0, 1 + shape * excess / scale)^
      ((-1/shape - 1) * a)
      c1 <- 1 / (scale^a * (1 + a + a * shape))
      c2 <- (1 + 1/a ) / scale^a
      div <- c1 - c2 * mean(y)
    }

    else
      div <- 1e6
    
    return(div)
  }
  
  opt <- optim(start, pddf, hessian = TRUE, ..., method = method)
  
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"
  
  shape <- opt$par[2]
  scale <- opt$par[1]
  
  param <- c(scale, shape)
  names(param) <- c("scale", "shape")
  
  std.err <- std.err.type <- var.cov <- corr <- NULL
  
  var.thresh <- FALSE
  
  list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = NULL, param = param,
       deviance = NULL, corr = corr, convergence = opt$convergence,
       counts = opt$counts, message = opt$message, threshold = threshold,
       nat = nat, pat = pat, data = x, exceed = exceed,
       scale = scale, var.thresh = var.thresh, est = "MDPD",
       opt.value = opt$value)
}



## This function comes from the evd package. The gpdmle function
## corresponds to the fpot function. Nevertheless, it was sligthly modified
## to simplify it. So, this function is a ligther version of fpot.
## So, I'm very gratefull to Alec Stephenson.

gpdmle <- function(x, threshold, start, ...,
                   std.err.type = "observed", corr = FALSE,
                   method = "BFGS", warn.inf = TRUE){

  if (all(c("observed", "expected", "none") != std.err.type))
    stop("``std.err.type'' must be one of 'observed', 'expected' or 'none'")
  
  nlpot <- function(scale, shape) { 
    -.C("gpdlik", exceed, nat, threshold, scale,
        shape, dns = double(1), PACKAGE = "simPop")$dns
  }
  
  nn <- length(x)
  
  threshold <- rep(threshold, length.out = nn)
  
  high <- (x > threshold) & !is.na(x)
  threshold <- as.double(threshold[high])
  exceed <- as.double(x[high])
  nat <- length(exceed)
  
  if(!nat) stop("no data above threshold")
  
  pat <- nat/nn
  param <- c("scale", "shape")
  
  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0)
    start$scale <- mean(exceed) - min(threshold)
    
    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nlpot)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlpot) <- c(f[m], f[-m])
  nllh <- function(p, ...) nlpot(p, ...)
  
  if(l > 1)
    body(nllh) <- parse(text = paste("nlpot(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if( warn.inf && do.call("nllh", start.arg) == 1e6 )
    warning("negative log-likelihood is infinite at starting values")
  
  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
    
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  if (std.err.type != "none"){
    
    tol <- .Machine$double.eps^0.5
    
    if(std.err.type == "observed") {
      
      var.cov <- qr(opt$hessian, tol = tol)
      if(var.cov$rank != ncol(var.cov$qr)){
        warning("observed information matrix is singular; passing std.err.type to ``expected''")
        obs.fish <- FALSE
        return
      }
      
      if (std.err.type == "observed"){
        var.cov <- try(solve(var.cov, tol = tol), silent = TRUE)

        if(!is.matrix(var.cov)){
          warning("observed information matrix is singular; passing std.err.type to ''none''")
          std.err.type <- "expected"
          return
        }

        else{
          std.err <- diag(var.cov)
          if(any(std.err <= 0)){
            warning("observed information matrix is singular; passing std.err.type to ``expected''")
            std.err.type <- "expected"
            return
          }
          
          std.err <- sqrt(std.err)
        
          if(corr) {
            .mat <- diag(1/std.err, nrow = length(std.err))
            corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr.mat) <- rep(1, length(std.err))
          }
          else {
            corr.mat <- NULL
          }
        }
      }
    }
    
    if (std.err.type == "expected"){
      
      shape <- opt$par[2]
      scale <- opt$par[1]
      a22 <- 2/((1+shape)*(1+2*shape))
      a12 <- 1/(scale*(1+shape)*(1+2*shape))
      a11 <- 1/((scale^2)*(1+2*shape))
      ##Expected Matix of Information of Fisher
      expFisher <- nat * matrix(c(a11,a12,a12,a22),nrow=2)

      expFisher <- qr(expFisher, tol = tol)
      var.cov <- solve(expFisher, tol = tol)
      std.err <- sqrt(diag(var.cov))
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      else
        corr.mat <- NULL
    }

    colnames(var.cov) <- nm
    rownames(var.cov) <- nm
    names(std.err) <- nm
  }

  else{
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }
  
  
  param <- c(opt$par, unlist(fixed.param))
  scale <- param["scale"]
  
  var.thresh <- !all(threshold == threshold[1])

  if (!var.thresh)
    threshold <- threshold[1]
  
  list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = unlist(fixed.param), param = param,
       deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
       counts = opt$counts, message = opt$message, threshold = threshold,
       nat = nat, pat = pat, data = x, exceed = exceed, scale = scale,
       var.thresh = var.thresh, est = "MLE", logLik = -opt$value,
       opt.value = opt$value, hessian = opt$hessian)
}

##Medians estimation for the GPD ( Peng, L. and Welsh, A. (2002) )
gpdmed <- function(x, threshold, start, tol = 10^-3, maxit = 500,
                   show.trace = FALSE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'med'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  nn <- length(x)
  
  threshold <- rep(threshold, length.out = nn)
  
  high <- (x > threshold) & !is.na(x)
  threshold <- as.double(threshold[high])
  exceed <- as.double(x[high])
  nat <- length(exceed)
  
  excess <- exceed - threshold
  
  if(!nat) stop("no data above threshold")
  
  pat <- nat/nn
  
  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0.1)
    start["scale"] <- mean(exceed) - min(threshold)
    
  }

  start <- c(scale = start$scale, shape = start$shape)
  iter <- 1
  
  trace <- round(start, 3)
  
  ##Definition of a function to solve
  f <- function(x, y){
    -log(x)/y - (1+y)/y^2 * (1 - x^y) + log(x + .5)/y +
      (1+y)/y^2 * (1 - (x+.5)^y)
  }
  
  
  while (iter < maxit){
    ##If we have a non feasible point, we move back to feasible region
    if ( (start[2] < 0) & (max(excess) >= (-start[1] / start[2])))
      start[2] <- -start[1] / max(excess) + .1
    
    r1 <- start[2] * median(excess) / (2^start[2] - 1) - start[1]
    
    a <- log( 1 + start[2] * excess / start[1] ) / start[2]^2
    b <- (1 + start[2]) * excess / (start[1]*start[2] +
                                    start[2]^2 * excess)
    
    if (start[2] <= -1)
      y1 <- .5
    
    else{
      opt <- uniroot(f, c(10^-12, .5), y = start[2])
      y1 <- opt$root
    }
    
    r2 <- median(a - b) + log(y1)/start[2] + (1 + start[2]) /
      start[2]^2 * (1 - y1^start[2])
    
    next.point <- c(r1, r2) + start
    
    if (sqrt(sum( (next.point - start)^2) ) < tol)
      break
    
    trace <- rbind(trace, next.point)
    iter <- iter + 1
    start <- next.point
    
  }
  
  if(iter == maxit) opt$convergence <- "iteration limit reached"
  
  else opt$convergence <- "successful"
  
  opt$counts <- iter - 1
  names(opt$counts) <- "function"
  
  shape <- start[2]
  scale <- start[1]
  
  param <- c(scale = scale, shape = shape)
  names(param) <- c("scale", "shape")
  
  std.err <- std.err.type <- var.cov <- corr <- NULL
  
  var.thresh <- FALSE
  
  
  if (show.trace){
    if (iter >= 2)
      rownames(trace) <- c("Init. Val.", 1:(iter-1))

    print(round(trace, 3))
  }
  
  list(fitted.values = param, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = NULL, param = param,
       deviance = NULL, corr = corr, convergence = opt$convergence,
       counts = opt$counts, message = opt$message, threshold = threshold,
       nat = nat, pat = pat, data = x, exceed = exceed,
       scale = scale, var.thresh = var.thresh, est = "MEDIANS",
       opt.value = opt$f.root)
  
}

