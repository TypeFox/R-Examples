## different leverage measures
## both currently have implementation that require
## n x n computations...hence these computations are done
## in the extractor functions (not the fitting function)

gleverage <- function(model, ...) {
  UseMethod("gleverage")
}

gleverage.betareg <- function(model, ...)
{
  ## extract response y and regressors X
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model, model = "mean") else model$x$mean
  z <- if(is.null(model$x)) model.matrix(model, model = "precision") else model$x$precision
  if(NCOL(x) < 1L) return(structure(rep.int(0, NROW(x)), .Names = rownames(x)))

  if(is.null(model$offset$mean)) model$offset$mean <- rep(0, NROW(x))
  if(is.null(model$offset$precision)) model$offset$precision <- rep(0, NROW(z))
  wts <- weights(model)
  if(is.null(wts)) wts <- 1
  ystar <- qlogis(y)
  
  ## extract coefficients
  beta <- model$coefficients$mean
  gamma <- model$coefficients$precision

  ## compute different types of "fitted" values
  eta <- as.vector(x %*% beta + model$offset[[1L]])
  phi_eta <- as.vector(z %*% gamma + model$offset[[2L]])
  mu <- model$link$mean$linkinv(eta)
  phi <- model$link$precision$linkinv(phi_eta)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ## compute first and second derivatives
  dmu <- as.vector(model$link$mean$mu.eta(eta))
  dphi <- as.vector(model$link$precision$mu.eta(phi_eta))

  d2mu <- switch(model$link$mean$name,
    "logit" = { dlogis(eta) * (1 - 2 * exp(eta)/(1 + exp(eta))) },
    "probit" = { -dnorm(eta) * eta },
    "cloglog" = { exp(-exp(eta)) * exp(eta) * (1 - exp(eta)) },
    "cauchit" = { -2 * eta/(1 + eta^2) * dcauchy(eta) },
    "log" = { mu },
    "loglog" = { exp(-exp(-eta)) * exp(-eta) * (exp(-eta) - 1) }
  )
  d2phi <- switch(model$link$precision$name,
    "identity" = 0,
    "log" = phi,
    "sqrt" = 2
  )
  
  ## compute L (below equation 5)
  a <- psi1 + psi2
  b <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
  v <- mu * (ystar - mustar) + digamma(phi) - digamma(phi * (1 - mu)) + log(1 - y)  
  Qbb <- phi * (phi * a * dmu^2 - (ystar - mustar) * d2mu)
  Qbg <- (phi * (mu * a - psi2) + mustar - ystar) * dmu * dphi
  Qgg <- b * dphi^2 - v * d2phi
  L <- crossprod(x,  Qbg * wts * z)
  L <- rbind(
    cbind(crossprod(x, Qbb * wts * x), L),
    cbind(t(L), crossprod(z, Qgg * wts * z))
  )

  ## compute D (below equation 11)
  D <- cbind(dmu * x, matrix(0, ncol = ncol(z), nrow = nrow(z)))

  ## compute Lty (below equation 11)
  Mb <- 1/(y * (1 - y))
  Mg <- (mu - y)/(y * (1 - y))
  Lty <- wts * t(cbind(phi * dmu * Mb * x, dphi * Mg * z))
  
  ## equation 11
  GL <- D %*% solve(L) %*% Lty
  return(wts * diag(GL))
}

hatvalues.betareg <- function(model, ...)
{
  ## extract response y and regressors X and Z
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model, model = "mean") else model$x$mean
  z <- if(is.null(model$x)) model.matrix(model, model = "precision") else model$x$precision
  if(NCOL(x) < 1L) return(structure(rep.int(0, NROW(x)), .Names = rownames(x)))
  
  if(is.null(model$offset$mean)) model$offset$mean <- rep(0, NROW(x))
  if(is.null(model$offset$precision)) model$offset$precision <- rep(0, NROW(z))
  wts <- weights(model)
  if(is.null(wts)) wts <- 1
  
  ## extract coefficients
  beta <- model$coefficients$mean
  gamma <- model$coefficients$precision

  ## compute different types of "fitted" values
  eta <- as.vector(x %*% beta + model$offset[[1L]])
  phi_eta <- as.vector(z %*% gamma + model$offset[[2L]])
  mu <- model$link$mean$linkinv(eta)
  phi <- model$link$precision$linkinv(phi_eta)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)

  ## compute w
  w <- wts * phi * (psi1 + psi2) * as.vector(model$link$mean$mu.eta(eta))^2

  ## compute (X'W X)^(-1)
  xwx1 <- chol2inv(qr.R(qr(sqrt(w) * x)))

  ## hat values
  ## (computations are bad: of order n x n)
  w * diag(x %*% xwx1 %*% t(x))
}
