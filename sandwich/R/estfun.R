estfun <- function(x, ...)
{
  UseMethod("estfun")
}

estfun.lm <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}

estfun.mlm <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  cf <- coef(x)
  rval <- lapply(1:NCOL(res), function(i) {
    rv <- as.vector(res[,i]) * wts * xmat
    colnames(rv) <- paste(colnames(cf)[i], colnames(rv), sep = ":")
    rv
  })  
  rval <- do.call("cbind", rval)
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(any(alias <- is.na(as.vector(cf)))) rval <- rval[, !alias, drop = FALSE]
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}

estfun.glm <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
    else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), na.rm = TRUE)
  rval <- wres * xmat / dispersion
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  res <- residuals(x, type = "pearson")
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}

estfun.rlm <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  psi <- function(z) x$psi(z) * z
  rval <- as.vector(psi(res/x$s)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}

estfun.polr <- function(x, ...)
{
  ## link processing
  mueta <- x$method
  if(mueta == "logistic") mueta <- "logit"
  mueta <- make.link(mueta)$mu.eta
  
  ## observations
  xmat <- model.matrix(x)[, -1L, drop = FALSE]
  n <- nrow(xmat)
  k <- ncol(xmat)
  m <- length(x$zeta)
  mf <- model.frame(x)
  y <- as.numeric(model.response(mf))
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, n)

  ## estimates  
  prob <- x$fitted.values[cbind(1:n, y)]
  xb <- if(k >= 1L) as.vector(xmat %*% x$coefficients) else rep(0, n)
  zeta <- x$zeta
  lp <- cbind(0, mueta(matrix(zeta, nrow = n, ncol = m, byrow = TRUE) - xb), 0)

  ## estimating functions
  rval <- matrix(0, nrow = n, ncol = k + m + 2L)
  if(k >= 1L) rval[, 1L:k] <- (-xmat * as.vector(lp[cbind(1:n, y + 1L)] - lp[cbind(1:n, y)]))
  rval[cbind(1:n, k + y)] <- -as.vector(lp[cbind(1:n, y)])
  rval[cbind(1:n, k + y + 1L)] <- as.vector(lp[cbind(1:n, y + 1L)])
  rval <- rval[, -c(k + 1L, k + m + 2L), drop = FALSE]
  rval <- w/prob * rval

  ## dimnames and return
  dimnames(rval) <- list(rownames(xmat), c(colnames(xmat), names(x$zeta)))
  return(rval)
}

estfun.clm <- function(x, ...)
{
  if(x$threshold != "flexible") stop("only flexible thresholds implemented at the moment")

  ## link processing
  mueta <- make.link(x$link)$mu.eta
  
  ## observations
  xmat <- model.matrix(x)
  if(length(xmat) > 1L) stop("estimating functions for scale regression not implemented yet")
  xmat <- xmat$X[, -1L, drop = FALSE]
  n <- nrow(xmat)
  k <- ncol(xmat)
  m <- length(x$alpha)
  mf <- model.frame(x)
  y <- as.numeric(model.response(mf))
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, n)

  ## estimates  
  prob <- x$fitted.values
  xb <- if(k >= 1L) as.vector(xmat %*% x$beta) else rep(0, n)
  zeta <- x$alpha
  lp <- cbind(0, mueta(matrix(zeta, nrow = n, ncol = m, byrow = TRUE) - xb), 0)

  ## estimating functions
  rval <- matrix(0, nrow = n, ncol = k + m + 2L)
  if(k >= 1L) rval[, 1L:k] <- (-xmat * as.vector(lp[cbind(1:n, y + 1L)] - lp[cbind(1:n, y)]))
  rval[cbind(1:n, k + y)] <- -as.vector(lp[cbind(1:n, y)])
  rval[cbind(1:n, k + y + 1L)] <- as.vector(lp[cbind(1:n, y + 1L)])
  rval <- rval[, -c(k + 1L, k + m + 2L), drop = FALSE]
  rval <- w/prob * rval

  ## dimnames, re-order and return
  dimnames(rval) <- list(rownames(xmat), c(colnames(xmat), names(x$alpha)))
  ix <- if(k >= 1L) c((k + 1L):(k + m), 1L:k) else 1L:m
  return(rval[, ix, drop = FALSE])
}

estfun.coxph <- function(x, ...)
{
  wts <- x$weights
  if(is.null(wts)) wts <- 1
  wts * residuals(x, type = "score", ...)
}

estfun.survreg <- function(x, ...)
{
  mf <- model.frame(x)
  xmat <- model.matrix(terms(x), mf)
  wts <- model.weights(mf)
  if(is.null(wts)) wts <- 1
  res <- residuals(x, type = "matrix")
  rval <- as.vector(res[,"dg"]) * wts * xmat
  if(NROW(x$var) > length(coef(x))) {
    rval <- cbind(rval, res[,"ds"])
    colnames(rval)[NCOL(rval)] <- "Log(scale)"
  }
  attr(rval, "assign") <- NULL
  
  return(rval)
}

estfun.nls <- function(x, ...)
{
  rval <- as.vector(x$m$resid()) * x$m$gradient()
  colnames(rval) <- names(coef(x))
  rval
}

estfun.hurdle <- function(x, ...) {
  ## extract data
  Y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  X <- model.matrix(x, model = "count")
  Z <- model.matrix(x, model = "zero")
  beta <- coef(x, model = "count")
  gamma <- coef(x, model = "zero")
  fulltheta <- x$theta

  offset <- x$offset
  if(is.list(offset)) {
    offsetx <- offset$count
    offsetz <- offset$zero
  } else {
    offsetx <- offset
    offsetz <- NULL
  }
  if(is.null(offsetx)) offsetx <- 0
  if(is.null(offsetz)) offsetz <- 0
  if(x$dist$zero == "binomial") linkobj <- make.link(x$link)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  Y1 <- Y > 0

  ## count component: working residuals
  eta <- as.vector(X %*% beta + offsetx)
  mu <- exp(eta)
  theta <- fulltheta["count"]

  wres_count <- as.numeric(Y > 0) * switch(x$dist$count,
    "poisson" = {
      (Y - mu) - exp(ppois(0, lambda = mu, log.p = TRUE) -
        ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)    
    },
    "geometric" = {
      (Y - mu * (Y + 1)/(mu + 1)) - exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
        pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) - log(mu + 1) + eta)
    },
    "negbin" = {
      (Y - mu * (Y + theta)/(mu + theta)) - exp(pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) +
	log(theta) - log(mu + theta) + eta)
    })
  
  ## zero component: working residuals
  eta <- as.vector(Z %*% gamma + offsetz)
  mu <- if(x$dist$zero == "binomial") linkobj$linkinv(eta) else exp(eta)
  theta <- fulltheta["zero"]

  wres_zero <- switch(x$dist$zero,
    "poisson" = {
      ifelse(Y1, exp(ppois(0, lambda = mu, log.p = TRUE) -
        ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta), -mu)
    },
    "geometric" = {
      ifelse(Y1, exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
        pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) - log(mu + 1) + eta), -mu/(mu + 1))
    },
    "negbin" = {
      ifelse(Y1, exp(pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) +
        log(theta) - log(mu + theta) + eta), -mu * theta/(mu + theta))
    },
    "binomial" = {
      ifelse(Y1, 1/mu, -1/(1-mu)) * linkobj$mu.eta(eta)
    })

  ## compute gradient from data
  rval <- cbind(wres_count * wts * X, wres_zero * wts * Z)
  colnames(rval) <- names(coef(x))
  rownames(rval) <- rownames(X)
  return(rval)
}

estfun.zeroinfl <- function(x, ...) {
  ## extract data
  Y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  X <- model.matrix(x, model = "count")
  Z <- model.matrix(x, model = "zero")
  beta <- coef(x, model = "count")
  gamma <- coef(x, model = "zero")
  theta <- x$theta

  offset <- x$offset
  if(is.list(offset)) {
    offsetx <- offset$count
    offsetz <- offset$zero
  } else {
    offsetx <- offset
    offsetz <- NULL
  }
  if(is.null(offsetx)) offsetx <- 0
  if(is.null(offsetz)) offsetz <- 0
  linkobj <- make.link(x$link)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  Y1 <- Y > 0

  eta <- as.vector(X %*% beta + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(Z %*% gamma + offsetz)
  muz <- linkobj$linkinv(etaz)

  ## density for y = 0
  clogdens0 <- switch(x$dist,
    "poisson" = -mu,
    "geometric" = dnbinom(0, size = 1, mu = mu, log = TRUE),
    "negbin" = dnbinom(0, size = theta, mu = mu, log = TRUE))
  dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

  ## working residuals  
  wres_count <- switch(x$dist,
    "poisson" = ifelse(Y1, Y - mu, -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(mu))),
    "geometric" = ifelse(Y1, Y - mu * (Y + 1)/(mu + 1), -exp(-log(dens0) +
      log(1 - muz) + clogdens0 - log(mu + 1) + log(mu))),
    "negbin" = ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), -exp(-log(dens0) +
      log(1 - muz) + clogdens0 + log(theta) - log(mu + theta) + log(mu))))
  wres_zero <- ifelse(Y1, -1/(1-muz) * linkobj$mu.eta(etaz),
    (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)

  ## compute gradient from data
  rval <- cbind(wres_count * wts * X, wres_zero * wts * Z)
  colnames(rval) <- names(coef(x))
  rownames(rval) <- rownames(X)
  return(rval)
}

estfun.mlogit <- function(x, ...)
{
  x$gradient
}
