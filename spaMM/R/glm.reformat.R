glm.reformat <- function (x, y, weights = rep(1, nobs), start, etastart, 
          mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
          control = list(), intercept = TRUE) 
{
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) 
    rownames(y)
  else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights)) 
    weights <- rep.int(1, nobs)
  if (is.null(offset)) 
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("'family' argument seems not to be a valid family object", 
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x)) 
    if.null
  else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  n <- NULL ## to avoid an R CMD check NOTE which cannot see that n will be set by eval(family$initialize)
  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta)) 
      stop("invalid linear predictor values in empty model", 
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu)) 
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  }
  else {
    coefold <- NULL
    eta <- if (!is.null(etastart)) 
      etastart
    else if (!is.null(start)) 
      if (length(start) != nvars) 
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                      nvars, paste(deparse(xnames), collapse = ", ")), 
             domain = NA)
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1L) 
        x * start
        else x %*% start)
    }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
      stop("cannot find valid starting values: please specify some", 
           call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    ###### replaces the main loop of glm.fit
    good <- weights > 0
    varmu <- variance(mu)[good]
    if (getRversion() > "3.1.0") {
      if (anyNA(varmu)) stop("NAs in V(mu)")
    } else {
      if (any(is.na(varmu))) stop("NAs in V(mu)")
    }
    if (any(varmu == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    #     if (all(!good)) {
    #       conv <- FALSE
    #       warning(gettextf("no observations informative at iteration %d", 
    #                        iter), domain = NA)
    #       break
    #     }
    z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
    w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
    #fit <- .Call(C_Cdqrls, x[good, , drop = FALSE] * 
    #               w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
    coef <- start
    names(coef) <- colnames(x)
    if (getRversion() > "3.1.0") {
      # .lm.fit is a wrapper for C_Cdqrls, but with check=TRUE
      fit <- .lm.fit(x[good, , drop = FALSE] * w, z * w, min(1e-07, control$epsilon/1000))
    } else {
      ##### reconstructing the C_Cdqrls call: 
      fit <- qr(x[good, , drop = FALSE] * w)
      fit$effects <- qr.qy(fit,z*w)
      #fit$coefficients <- qr.solve(fit,z*w) # not useful: not used below nor in return value
    }
    coef[fit$pivot] <- unique(etastart) ## keep names
    ######
    eta <- etastart
    mu <- linkinv(eta)
    dev <- sum(dev.resids(y, mu, weights))
    if (fit$rank < nvars) 
      coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    residuals <- (y - mu)/mu.eta(eta)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    xxnames <- xnames[fit$pivot]
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY) 
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
                                                                sum(good) - fit$rank))
  wtdmu <- if (intercept) 
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) 
    0
  else fit$rank
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
  # derogative comment from glm.fit source:
  ##     ^^ is only initialize()d for "binomial" [yuck!]
  # n is in the envir of binomial() and set by binomial()$initialize() called above by eval(family$initialize) 
  list(coefficients = coef, residuals = residuals, fitted.values = mu, 
       effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
       rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
                                                     "qraux", "pivot", "tol")], class = "qr"), family = family, 
       linear.predictors = eta, deviance = dev, aic = aic.model, 
       null.deviance = nulldev, iter = 1L, weights = wt, prior.weights = weights, 
       df.residual = resdf, df.null = nulldf, y = y, converged = NULL, 
       boundary = NULL)
}