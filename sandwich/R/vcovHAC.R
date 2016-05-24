## vcovHAC() is the general workhorse for HAC estimation
## although the essential part is now moved to meatHAC()
## interfacing the sandwich() function

vcovHAC <- function(x, ...) {
  UseMethod("vcovHAC")
}

vcovHAC.default <- function(x, order.by = NULL, prewhite = FALSE,
  weights = weightsAndrews, adjust = TRUE, diagnostics = FALSE,
  sandwich = TRUE, ar.method = "ols", data = list(), ...)
{
  rval <- meatHAC(x, order.by = order.by, prewhite = prewhite,
                  weights = weights, adjust = adjust, diagnostics = diagnostics,
		  ar.method = ar.method, data = data)
				  
  if(sandwich) {
    diagn <- attr(rval, "diagnostics")
    rval <- sandwich(x, meat = rval, ...)
    attr(rval, "diagnostics") <- diagn
  }

  return(rval)
}

meatHAC <- function(x, order.by = NULL, prewhite = FALSE,
  weights = weightsAndrews, adjust = TRUE, diagnostics = FALSE,
  ar.method = "ols", data = list())
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  prewhite <- as.integer(prewhite)

  umat <- estfun(x)[, , drop = FALSE]
  if(is.zoo(umat)) umat <- as.matrix(coredata(umat))
  n.orig <- n <- nrow(umat)
  k <- ncol(umat)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  umat <- umat[index, , drop = FALSE]

  if(prewhite > 0) {
    var.fit <- try(ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method))
    if(inherits(var.fit, "try-error")) stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    if(k > 1) D <- solve(diag(ncol(umat)) - apply(var.fit$ar, 2:3, sum))
      else D <- as.matrix(1/(1 - sum(var.fit$ar)))
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }

  if(is.function(weights))
    weights <- weights(x, order.by = order.by, prewhite = prewhite, ar.method = ar.method, data = data)

  if(length(weights) > n) {
    warning("more weights than observations, only first n used")
    weights <- weights[1:n]
  }
 
  utu <- 0.5 * crossprod(umat) * weights[1]
  wsum<-n*weights[1]/2
  w2sum<-n*weights[1]^2/2

  if(length(weights) > 1) {
    for (ii in 2:length(weights)) {
      utu <- utu + weights[ii] * crossprod(umat[1:(n-ii+1),,drop=FALSE], umat[ii:n,,drop=FALSE])
      wsum <- wsum + (n-ii+1) * weights[ii]
      w2sum <- w2sum + (n-ii+1) * weights[ii]^2
    }
  }

  utu <- utu + t(utu)

  ## Andrews: multiply with df n/(n-k)
  if(adjust) utu <- n.orig/(n.orig-k) * utu
  
  if(prewhite > 0) {
    utu <- crossprod(t(D), utu) %*% t(D)
  }
  
  wsum <- 2*wsum
  w2sum <- 2*w2sum
  bc <- n^2/(n^2 - wsum)
  df <- n^2/w2sum

  rval <- utu/n.orig

  if(diagnostics)  attr(rval, "diagnostics") <- list(bias.correction = bc, df = df)
  return(rval)
}


## weightsAndrews() and bwAndrews() implement the HAC estimation
## procedure described in Andrews (1991) and Andrews & Monahan (1992)
## kernHAC() is the convenience interface.
## (Note, that bwNeweyWest() can also be used with weightsAndrews())

weightsAndrews <- function(x, order.by = NULL, bw = bwAndrews,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  prewhite = 1, ar.method = "ols", tol = 1e-7, data = list(), verbose = FALSE, ...)
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  kernel <- match.arg(kernel)
  if(is.function(bw))
    bw <- bw(x, order.by = order.by, kernel = kernel,
      prewhite = prewhite, data = data, ar.method = ar.method, ...)
  if(verbose) cat(paste("\nBandwidth chosen:", format(bw), "\n"))
      
  n <- NROW(estfun(x)) - as.integer(prewhite)
  
  weights <- kweights(0:(n-1)/bw, kernel = kernel)
  weights <- weights[1:max(which(abs(weights) > tol))]
  return(weights)
}

bwAndrews <- function(x, order.by = NULL, kernel = c("Quadratic Spectral", "Truncated",
  "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
  weights = NULL, prewhite = 1, ar.method = "ols", data = list(), ...)
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  kernel <- match.arg(kernel)
  approx <- match.arg(approx)
  prewhite <- as.integer(prewhite)

  umat <- estfun(x)[,, drop = FALSE]
  if(is.zoo(umat)) umat <- as.matrix(coredata(umat))
  n <- nrow(umat)
  k <- ncol(umat)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }

  umat <- umat[index, , drop = FALSE]

  ## compute weights (try to set the intercept weight to 0)
  #### could be ignored by using: weights = 1
  
  if(is.null(weights)) {
    weights <- rep(1, k)
    unames <- colnames(umat)
    if(!is.null(unames) && "(Intercept)" %in% unames)
      weights[which(unames == "(Intercept)")] <- 0
    else {
      res <- try(as.vector(rowMeans(estfun(x)/model.matrix(x), na.rm = TRUE)), silent = TRUE)
      if(inherits(res, "try-error")) res <- try(residuals(x), silent = TRUE)
      if(!inherits(res, "try-error")) weights[which(colSums((umat - res)^2) < 1e-16)] <- 0
    }
    if(isTRUE(all.equal(weights, rep(0, k)))) weights <- rep(1, k)
  } else {
    weights <- rep(weights, length.out = k)
  }
  if(length(weights) < 2) weights <- 1

  if(prewhite > 0) {
    var.fit <- ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    if(inherits(var.fit, "try-error")) stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite ##??
  }

  if(approx == "AR(1)") {
    fitAR1 <- function(x) {
      rval <-  ar(x, order.max = 1, aic = FALSE, method = "ols")
      rval <- c(rval$ar, sqrt(rval$var.pred))
      names(rval) <- c("rho", "sigma")
      return(rval)
    }

    ar.coef <- apply(umat, 2, fitAR1)

    denum <- sum(weights * (ar.coef["sigma",]/(1-ar.coef["rho",]))^4)
    alpha2 <- sum(weights * 4 * ar.coef["rho",]^2 * ar.coef["sigma",]^4/(1-ar.coef["rho",])^8) / denum
    alpha1 <- sum(weights * 4 * ar.coef["rho",]^2 * ar.coef["sigma",]^4/((1-ar.coef["rho",])^6 * (1+ar.coef["rho",])^2)) / denum

  } else {

    fitARMA11 <- function(x) {
      rval <-  arima(x, order = c(1, 0, 1), include.mean = FALSE)
      rval <- c(rval$coef, sqrt(rval$sigma2))
      names(rval) <- c("rho", "psi", "sigma")
      return(rval)
    }

    arma.coef <- apply(umat, 2, fitARMA11)

    denum <- sum(weights * ((1 + arma.coef["psi",]) * arma.coef["sigma",]/(1-arma.coef["rho",]))^4)
    alpha2 <- sum(weights * 4 * ((1 + arma.coef["rho",] * arma.coef["psi",]) * (
                                  arma.coef["rho",] + arma.coef["psi",]))^2 * arma.coef["sigma",]^4/
				 (1-arma.coef["rho",])^8) / denum
    alpha1 <- sum(weights * 4 * ((1 + arma.coef["rho",] * arma.coef["psi",]) * (
                                  arma.coef["rho",] + arma.coef["psi",]))^2 * arma.coef["sigma",]^4/
                                 ((1-arma.coef["rho",])^6 * (1+arma.coef["rho",])^2)) / denum
  }

  rval <- switch(kernel,
    "Truncated"          = {0.6611 * (n * alpha2)^(1/5)},
    "Bartlett"           = {1.1447 * (n * alpha1)^(1/3)},
    "Parzen"             = {2.6614 * (n * alpha2)^(1/5)},
    "Tukey-Hanning"      = {1.7462 * (n * alpha2)^(1/5)},
    "Quadratic Spectral" = {1.3221 * (n * alpha2)^(1/5)})

  return(rval)  
}

kernHAC <- function(x, order.by = NULL, prewhite = 1, bw = bwAndrews,
  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
  approx = c("AR(1)", "ARMA(1,1)"), adjust = TRUE, diagnostics = FALSE, sandwich = TRUE,
  ar.method = "ols", tol = 1e-7, data = list(), verbose = FALSE, ...)
{
  myweights <- function(x, order.by = NULL, prewhite = FALSE, ar.method = "ols", data = list())
    weightsAndrews(x, order.by = order.by, prewhite = prewhite, bw = bw,
                   kernel = kernel, approx = approx, ar.method = ar.method, tol = tol,
		   data = data, verbose = verbose, ...)
  vcovHAC(x, order.by = order.by, prewhite = prewhite,
    weights = myweights, adjust = adjust, diagnostics = diagnostics,
    sandwich = sandwich, ar.method = ar.method, data = data)
}



## weightsLumley() implements the WEAVE estimators from 
## Lumley & Heagerty (1999)
## weave() is a convenience interface

weightsLumley <- function(x, order.by = NULL, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, tol = 1e-7, data = list(), ...)
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  method <- match.arg(method)
  res <- residuals(x, "response") #FIXME# available for which models?
  n <- length(res)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  res <- res[index]

  rhohat <- acf(res)

  switch(method,
  "truncate" = {
    if(is.null(C)) C <- 4
    lag <- max((1:length(rhohat))[rhohat^2*n > C])
    weights <- rep(1, lag)
  },
  "smooth" = {
    if(is.null(C)) C <- 1
    weights <- C * n * rhohat^2
    weights <- ifelse(weights > 1, 1, weights)
    weights <- weights[1:max(which(abs(weights) > tol))]
  })
  
  return(weights)
}



weave <- function(x, order.by = NULL, prewhite = FALSE, C = NULL,
  method = c("truncate", "smooth"), acf = isoacf, adjust = FALSE,
  diagnostics = FALSE, sandwich = TRUE, tol = 1e-7, data = list(), ...)
{
  myweights <- function(x, order.by = NULL, prewhite = FALSE, data = list(), ...)
    weightsLumley(x, order.by = order.by, prewhite = prewhite, C = C,
                   method = method, acf = acf, tol = tol, data = data)
  vcovHAC(x, order.by = order.by, prewhite = prewhite,
    weights = myweights, adjust = adjust, diagnostics = diagnostics,
    sandwich = sandwich, data = data)
}


## bwNeweyWest() implements the procedure from Newey & West (1994)
## It works for Bartlett/Parzen/QS kernels and can thus be passed
## to weightsAndrews() and kernHAC() respectively.
## A convenience interface NeweyWest() to only the Bartlett kernel
## is also available.

bwNeweyWest <- function(x, order.by = NULL, kernel = c("Bartlett", "Parzen",
  "Quadratic Spectral", "Truncated", "Tukey-Hanning"), weights = NULL, prewhite = 1,
  ar.method = "ols", data = list(), ...)
{
  ## ensure that NAs are omitted
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"

  kernel <- match.arg(kernel)
  if(kernel %in% c("Truncated", "Tukey-Hanning"))
    stop(paste("Automatic bandwidth selection only available for ", 
      dQuote("Bartlett"), ", ", dQuote("Parzen"), " and ", dQuote("Quadratic Spectral"),
      " kernel. Use ", sQuote("bwAndrews"), " instead.", sep = ""))
  prewhite <- as.integer(prewhite)

  umat <- estfun(x)[,, drop = FALSE]
  if(is.zoo(umat)) umat <- as.matrix(coredata(umat))
  n <- nrow(umat)
  k <- ncol(umat)

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }

  umat <- umat[index, , drop = FALSE]

  ## compute weights (try to set the intercept weight to 0)
  #### could be ignored by using: weights = 1
  
  if(is.null(weights)) {
    weights <- rep(1, k)
    unames <- colnames(umat)
    if(!is.null(unames) && "(Intercept)" %in% unames)
      weights[which(unames == "(Intercept)")] <- 0
    else {
      res <- try(as.vector(rowMeans(estfun(x)/model.matrix(x), na.rm = TRUE)), silent = TRUE)
      if(inherits(res, "try-error")) res <- try(residuals(x), silent = TRUE)
      if(!inherits(res, "try-error")) weights[which(colSums((umat - res)^2) < 1e-16)] <- 0      
    }
    if(all(weights <= 0)) weights <- rep(1, length.out = k)
  } else {
    weights <- rep(weights, length.out = k)
  }
  if(length(weights) < 2) weights <- 1

  ## select lag truncation according to Table II C. from Newey & West (1994)
  mrate <- switch(kernel, 
    "Bartlett"           = 2/9,
    "Parzen"             = 4/25,
    "Quadratic Spectral" = 2/25)
  m <- floor(ifelse(prewhite > 0, 3, 4) * (n/100)^mrate)

  if(prewhite > 0) {
    var.fit <- ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    if(inherits(var.fit, "try-error")) stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }

  ## compute weighted variances
  hw <- umat %*% weights
  sigmaj <- function(j) sum(hw[1:(n-j)] * hw[(j+1):n])/n
  sigma <- sapply(0:m, sigmaj)
  s0 <- sigma[1] + 2*sum(sigma[-1])
  s1 <- 2 * sum(1:m * sigma[-1])
  s2 <- 2 * sum((1:m)^2 * sigma[-1])
  
  ## use parameters as in Table I B.
  ## choose 1/(2*q + 1)
  qrate <- 1/(2 * ifelse(kernel == "Bartlett", 1, 2) + 1)
  ## compute gamma
  rval <- switch(kernel,
    "Bartlett"           = { 1.1447 * ((s1/s0)^2)^qrate },
    "Parzen"             = { 2.6614 * ((s2/s0)^2)^qrate },
    "Quadratic Spectral" = { 1.3221 * ((s2/s0)^2)^qrate })
  ## compute bandwidth
  rval <- rval * (n + prewhite)^qrate

  ## rval is not truncated. This is done in NeweyWest(),
  ## but bwNeweyWest() can also be used without truncation.
  
  return(rval)  
}

NeweyWest <- function(x, lag = NULL,
  order.by = NULL, prewhite = TRUE, adjust = FALSE, 
  diagnostics = FALSE, sandwich = TRUE, ar.method = "ols", data = list(),
  verbose = FALSE)
{
  if(is.null(lag)) lag <- floor(bwNeweyWest(x, 
    order.by = order.by, prewhite = prewhite,
    ar.method = ar.method, data = data))
  if(verbose) cat(paste("\nLag truncation parameter chosen:", lag, "\n"))
  
  myweights <- seq(1, 0, by = -(1/(lag + 1)))
  vcovHAC(x, order.by = order.by, prewhite = prewhite,
    weights = myweights, adjust = adjust, diagnostics = diagnostics,
    sandwich = sandwich, ar.method = ar.method, data = data)
}

