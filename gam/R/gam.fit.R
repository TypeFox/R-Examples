"gam.fit" <-
  function (x, y, smooth.frame, weights = rep(1, nobs), start = NULL, 
            etastart = NULL, mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = gam.control()) 
{
  ynames <- if (is.matrix(y)) 
    dimnames(y)[[1]]
  else names(y)
  xnames <- dimnames(x)[[2]]
  nobs <- NROW(y)
  nvars <- ncol(x)
  maxit <- control$maxit
  bf.maxit <- control$bf.maxit
  epsilon <- control$epsilon
  bf.epsilon <- control$bf.epsilon
  trace <- control$trace
  digits <- -log10(epsilon) + 1
  if (is.null(weights)) 
    weights <- rep.int(1, nobs)
  if (is.null(offset)) 
    offset <- rep.int(0, nobs)
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
    validmu <- function(mu) TRUE
  eval(family$initialize)
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  eta <- if (!is.null(etastart)) 
    etastart
  else if (!is.null(start)) 
    if (length(start) != nvars) 
      stop("Length of start should equal ", nvars, " and correspond to initial coefs for ", 
           deparse(xnames))
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1) 
                         x * start
      else x %*% start)
    }
  else family$linkfun(mustart)
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) 
    stop("Can't find valid starting values: please specify some")
  new.dev <- sum(dev.resids(y, mu, weights))
  a <- attributes(attr(smooth.frame, "terms"))
  smoothers <- a$specials
  if (length(smoothers) > 0) {
    smoothers <- smoothers[sapply(smoothers, length) > 0]
    for (i in seq(along = smoothers)) {
      tt <- smoothers[[i]]
      ff <- apply(a$factors[tt, , drop = FALSE], 2, any)
      smoothers[[i]] <- if (any(ff)) 
        seq(along = ff)[a$order == 1 & ff]
      else NULL
    }
  }
  if (length(smoothers) > 0) {
    smooth.labels <- a$term.labels[unlist(smoothers)]
    assignx <- attr(x, "assign")
    assignx <- assign.list(assignx, a$term.labels)
    which <- assignx[smooth.labels]
    if (length(smoothers) > 1) 
      bf <- "general.wam"
    else {
      sbf <- match(names(smoothers), gam.wlist, FALSE)
      bf <- if (sbf) 
        paste(gam.wlist[sbf], "wam", sep = ".")
      else "general.wam"
    }
    bf.call <- parse(text = paste(bf, "(x, z, wz, fit$smooth, which, fit$smooth.frame,bf.maxit,bf.epsilon, trace)", 
                       sep = ""))[[1]]
    s <- matrix(0, length(y), length(which))
    dimnames(s) <- list(names(y), names(which))
    fit <- list(smooth = s, smooth.frame = smooth.frame)
  }
  else {
    bf.call <- expression(lm.wfit(x, z, wz, method = "qr", 
        singular.ok = TRUE))
    bf <- "lm.wfit"
  }
  old.dev <- 10 * new.dev + 10
  for (iter in 1:maxit) {
    good <- weights > 0
    varmu <- variance(mu)
    if (any(is.na(varmu[good]))) 
      stop("NAs in V(mu)")
    if (any(varmu[good] == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    z <- eta - offset
    z[good] <- z[good] + (y - mu)[good]/mu.eta.val[good]
    wz <- weights
    wz[!good] <- 0
    wz[good] <- wz[good] * mu.eta.val[good]^2/varmu[good]
    fit <- eval(bf.call)
    eta <- fit$fitted.values + offset
    mu <- linkinv(eta)
    old.dev <- new.dev
    new.dev <- sum(dev.resids(y, mu, weights))
    if (trace) 
      cat("GAM ", bf, " loop ", iter, ": deviance = ", 
          format(round(new.dev, digits)), " \n", sep = "")
    if (is.na(new.dev)) {
      one.more <- FALSE
      warning("iterations terminated prematurely because of singularities")
    }
    else one.more <- abs(old.dev - new.dev)/(old.dev + 0.1) > 
      epsilon
    if (!one.more) 
      break
  }
  fitqr <- fit$qr
  xxnames <- xnames[fitqr$pivot]
  nr <- min(sum(good), nvars)
  if (nr < nvars) {
    Rmat <- diag(nvars)
    Rmat[1:nr, 1:nvars] <- fitqr$qr[1:nr, 1:nvars]
  }
  else Rmat <- fitqr$qr[1:nvars, 1:nvars]
  Rmat <- as.matrix(Rmat)
  Rmat[row(Rmat) > col(Rmat)] <- 0
  dimnames(Rmat) <- list(xxnames, xxnames)
  names(fit$residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  fit$additive.predictors <- eta
  fit$fitted.values <- mu
  names(fit$weights) <- ynames
  names(fit$effects) <- c(xxnames[seq(len = fitqr$rank)], rep.int("", 
                                        sum(good) - fitqr$rank))
  if (length(fit$smooth) > 0) 
    fit$smooth.frame <- smooth.frame[smooth.labels]
  wtdmu <- if (a$intercept) 
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(a$intercept)
  rank <- n.ok - fit$df.residual
  aic.model <- aic(y, nobs, mu, weights, new.dev) + 2 * rank
  if (!is.null(fit$smooth)) {
    nonzeroWt <- (wz > 0)
    nl.chisq <-  gam.nlchisq(fit$qr, fit$residuals, wz, fit$smooth)
  }
  else nl.chisq <- NULL
  fit <- c(fit, list(R = Rmat, rank = fitqr$rank, family = family, 
                     deviance = new.dev, aic = aic.model, null.deviance = nulldev, 
                     iter = iter, prior.weights = weights, y = y, df.null = nulldf, 
                     nl.chisq = nl.chisq))
  fit
}

