crch <- function(formula, data, subset, na.action, weights, offset, 
  link.scale = c("log", "identity", "quadratic"), 
  dist = c("gaussian", "logistic", "student"), df = NULL,
  left = -Inf, right = Inf, truncated = FALSE , control = crch.control(...),
  model = TRUE, x = FALSE, y = FALSE, ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data, dot = control$dot)
  mtX <- terms(formula, data = data, rhs = 1L, dot = control$dot)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L, dot = control$dot))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## obtain correct subset of predvars/dataClasses to terms
  .add_predvars_and_dataClasses <- function(terms, model.frame) {
    ## original terms
    rval <- terms
    ## terms from model.frame
    nval <- if(inherits(model.frame, "terms")) model.frame else terms(model.frame, dot = control$dot)

    ## associated variable labels
    ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1]
    nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1]
    if(!all(ovar %in% nvar)) stop(
      paste("The following terms variables are not part of the model.frame:",
      paste(ovar[!(ovar %in% nvar)], collapse = ", ")))
    ix <- match(ovar, nvar)
  
    ## subset predvars
    if(!is.null(attr(rval, "predvars"))) 
      warning("terms already had 'predvars' attribute, now replaced")
    attr(rval, "predvars") <- attr(nval, "predvars")[1L + c(0L, ix)]

    ## subset dataClasses
    if(!is.null(attr(rval, "dataClasses"))) 
      warning("terms already had 'dataClasses' attribute, now replaced")
    attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]
  
    return(rval)
  }
  mt  <- .add_predvars_and_dataClasses(mt,  mf)
  mtX <- .add_predvars_and_dataClasses(mtX, mf)
  mtZ <- .add_predvars_and_dataClasses(mtZ, mf)

  ## link
  if(is.character(link.scale)) link.scale <- match.arg(link.scale)


  ## distribution
  if(is.character(dist)) dist <- match.arg(dist)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(identical(dist, "student")) {
    if(!is.null(df) && df <= 0) stop("'df' must be positive")
    if(!is.null(df) && !is.finite(df)) dist <- "gaussian"
  }

  ## convenience variables
  n <- length(Y)


  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## in location part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## in scale part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for location)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(location = offsetX, scale = offsetZ)
 

  ## call the actual workhorse: crch.fit() or crch.boost()
  fit <- control$fit
  control$fit <- NULL
  rval <- do.call(fit, list(x = X, y = Y, z = Z, left = left, right = right, 
      link.scale = link.scale, dist = dist, df = df, weights = weights, 
      offset = offset, control = control, truncated = truncated))
  

  ## further model information
  rval$call <- if(length(control$call)) control$call else cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), 
    scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z)

  return(rval)
}


trch <- function(formula, data, subset, na.action, weights, offset,
  link.scale = c("log", "identity", "quadratic"), 
  dist = c("gaussian", "logistic", "student"), df = NULL,
  left = -Inf, right = Inf, truncated = TRUE , control = crch.control(...),
  model = TRUE, x = FALSE, y = FALSE, ...) 
{
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- as.name("crch")
  cl2$truncated <- truncated
  cl2$call <- as.name("cl")
  
  eval(cl2)
}

crch.control <- function(method = "BFGS", maxit = NULL, 
  hessian = NULL, trace = FALSE, start = NULL, dot = "separate", ...)
{
  if(method == "boosting") {
    if(is.null(maxit)) maxit <- 100
    rval <- crch.boost(dot = dot, start = start, maxit = maxit, ...)
  } else {
    if(is.null(maxit)) maxit <- 5000
    rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, 
      start = start, dot = dot, fit = "crch.fit")
    rval <- c(rval, list(...))
    if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
    rval$fnscale <- 1
    if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  }
  rval
}


crch.fit <- function(x, z, y, left, right, truncated = FALSE, 
  dist = "gaussian", df = NULL, link.scale = "log",
  weights = NULL, offset = NULL, control = crch.control()) 
{
  ## response and regressor matrix
  n <- NROW(x)  
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  dfest <- identical(dist, "student") & is.null(df)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(location = offset, scale = rep.int(0, n))
  if(is.null(z)) {
    q <- 1L
    z <- matrix(1, ncol = q, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
  } else {
    q <- NCOL(z)
    if(q < 1L) stop("scale regression needs to have at least one parameter")
  }

  ## control parameters
  ocontrol <- control
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  control$method <- control$hessian <- control$start <- NULL

  

  if(is.character(dist)){
    ## distribution functions
    if(truncated) {
      ddist2 <- switch(dist, 
        "student"  = dtt, "gaussian" = dtnorm, "logistic" = dtlogis)
      sdist2 <- switch(dist, 
        "student"  = stt, "gaussian" = stnorm, "logistic" = stlogis)
      hdist2 <- switch(dist, 
        "student"  = htt, "gaussian" = htnorm, "logistic" = htlogis)
    } else {
      ddist2 <- switch(dist, 
        "student"  = dct, "gaussian" = dcnorm, "logistic" = dclogis)
      sdist2 <- switch(dist, 
        "student"  = sct, "gaussian" = scnorm, "logistic" = sclogis)
      hdist2 <- switch(dist, 
        "student"  = hct, "gaussian" = hcnorm, "logistic" = hclogis)
    }
    ddist <- if(dist == "student") ddist2 else function(..., df) ddist2(...)
    sdist <- if(dist == "student") sdist2 else function(..., df) sdist2(...)
    hdist <- if(dist == "student") hdist2 else function(..., df) hdist2(...)


  } else { 
    ## for user defined distribution (requires list with ddist, sdist (optional)
    ## and hdist (optional), ddist, sdist, and hdist must be functions with
    ## arguments x, mean, sd, df, left, right, and log)
    ddist <- dist$ddist
    sdist <- if(is.null(dist$sdist)) NULL else  dist$sdist
    if(is.null(dist$hdist)) {
      if(hessian == FALSE) warning("no analytic hessian available. Hessian is set to TRUE and numerical Hessian from optim is employed")
      hessian <- TRUE     
    } else hdist <- dist$hdist 
    dist <- "user defined"
  }

  ## analytic or numeric Hessian
  if(is.null(hessian)) {
    hessian <- dfest
    returnvcov <- TRUE  # vcov is not computed when hessian == FALSE
  } else {
    returnvcov <- hessian
  } 

  ## link
  if(is.character(link.scale)) {
    linkstr <- link.scale
    if(linkstr != "quadratic") {
      linkobj <- make.link(linkstr)
      linkobj$dmu.deta <- switch(linkstr, 
        "identity" = function(eta) rep.int(0, length(eta)), 
        "log" = function(eta) pmax(exp(eta), .Machine$double.eps))
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) mu^2,
        linkinv = function(eta) sqrt(eta),
        mu.eta = function(eta) 1/2/sqrt(eta),
        dmu.deta = function(eta) -1/4/sqrt(eta^3),
        valideta = function(eta) TRUE,
        name = "quadratic"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link.scale
    linkstr <- link.scale$name
    if(is.null(linkobj$dmu.deta) & !hessian) {
      warning("link.scale needs to provide dmu.deta component for analytical Hessian. Numerical Hessian is employed.")
      hessian <- TRUE
    }
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta


  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, y, w = weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    gamma <- c(linkfun(sqrt(sum(weights * auxreg$residuals^2)/
      auxreg$df.residual)), rep(0, ncol(z) - 1))
    start <- if(dfest) c(beta, gamma, log(10)) else c(beta, gamma)
  }
  if(is.list(start)) start <- do.call("c", start) 
  if(length(start) > k + q + dfest) {
    warning(paste("too many entries in start! only first", k + q + dfest, "entries are considered"))
    start <- start[1: (k + q + dfest)]
  }
  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = q) + k]
    delta <- if(dfest) tail(par, 1) else NULL
    mu <- drop(x %*% beta) + offset[[1L]]
    zgamma <- drop(z %*% gamma) + offset[[2L]]
    sigma <- linkinv(zgamma)
    df <- if(dfest) exp(delta) else df
    list(
      beta = beta,
      gamma = gamma,
      delta = delta,
      mu = mu,
      zgamma = zgamma,
      sigma = sigma,
      df = df
    )
  }
  ## objective function
  loglikfun <- function(par) {
    fit <- fitfun(par)
    ll <- with(fit,  
        ddist(y, mu, sigma, df = df, left = left, right = right, log = TRUE))
    if(any(!is.finite(ll))) NaN else -sum(weights * ll)  
  }
  ## functions to evaluate gradients and hessian
  if(dfest | is.null(sdist)) {
    gradfun <- NULL
  } else { 
    gradfun <- function(par, type = "gradient") {
      fit <- fitfun(par)
      grad <- with(fit, 
        sdist(y, mu, sigma, df = df, left = left, right = right))
      grad <- cbind(grad[,1]*x, grad[,2] * mu.eta(fit$zgamma) * z)
      return(-colSums(weights * grad))
    }
    hessfun <- function(par) {
      fit <- fitfun(par)
      hess <- with(fit, hdist(y, mu, sigma, left = left, right = right,
        df = df, which = c("mu", "sigma", "mu.sigma", "sigma.mu")))
      grad <- with(fit, sdist(y, mu, sigma, left = left, right = right, 
        df = df, which = "sigma"))
      hess[, "d2sigma"] <- hess[, "d2sigma"]*mu.eta(fit$zgamma)^2 + grad*dmu.deta(fit$zgamma)
      hess[, "dmu.dsigma"] <- hess[, "dsigma.dmu"] <- hess[, "dmu.dsigma"]*mu.eta(fit$zgamma)
      hess <- weights*hess
      hessmu <- crossprod(hess[,"d2mu"]*x, x)
      hessmusigma <- crossprod(hess[,"dmu.dsigma"]*x, z)
      hesssigmamu <- crossprod(hess[,"dsigma.dmu"]*z, x)
      hesssigma <- crossprod(hess[,"d2sigma"]*z, z)
      -cbind(rbind(hessmu, hesssigmamu), rbind(hessmusigma, hesssigma))
    }
  }
  opt <- suppressWarnings(optim(par = start, fn = loglikfun, gr = gradfun,
    method = method, hessian = hessian, control = control))
  if(opt$convergence > 0) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }
  par <- opt$par
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  delta <- fit$delta
  mu <- fit$mu
  sigma <- fit$sigma
  vcov <- if(returnvcov) {
    if (hessian) solve(as.matrix(opt$hessian)) 
    else solve(hessfun(par))
  } else matrix(NA, k+q+dfest, n+k+dfest)
  ll <- -opt$value
  df <- if(dfest) exp(delta) else df

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  if (returnvcov) {
    colnames(vcov) <- rownames(vcov) <- c(
      colnames(x),
      paste("(scale)", colnames(z), sep = "_"),
      if(dfest) "(Log(df))")
  }

  rval <- list(
    coefficients = list(location = beta, scale = gamma, df = delta),
    df = df,
    residuals = y - mu,
    fitted.values = list(location = mu, scale = sigma),
    dist = dist,
    cens = list(left = left, right = right),
    optim = opt,  
    method = method,
    control = ocontrol,
    start = start,  
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(location = if(identical(offset[[1L]], rep.int(0, n))) NULL else 
      offset[[1L]],
    scale = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    loglik = ll,
    vcov = vcov,
    link = list(scale = linkobj),
    truncated = truncated,
    converged = converged,
    iterations = as.vector(tail(na.omit(opt$counts), 1))
  )
  
  class(rval) <- "crch"
  return(rval)
}

print.crch <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients$location)) {
      cat(paste("Coefficients (location model):\n", sep = ""))
       print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients (in location model)\n\n")
    if(length(x$coefficients$scale)) {
      cat(paste("Coefficients (scale model with ", x$link$scale$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat(paste("Distribution: ", x$dist, "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }

  invisible(x)
}


summary.crch <- function(object, ...)
{
  ## residuals
  object$residuals <- object$residuals/object$fitted.values$scale

  ## extend coefficient table
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if(length(object$coefficients$df)) {
    cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE], 
      df = cf[nrow(cf), , drop = FALSE])
    rownames(cf$df) <- names(object$coefficients$df)
  } else {
    cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE])
  }
  rownames(cf$location) <- names(object$coefficients$location)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.crch"
  object
}


print.summary.crch <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$location)) {
      cat(paste("\nCoefficients (location model):\n", sep = ""))
      printCoefmat(x$coefficients$location, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in location model)\n")

    if(NROW(x$coefficients$scale)) {
      cat(paste("\nCoefficients (scale model with ", x$link$scale$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients ( in scale model)\n")

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat(paste("\nDistribution: ", x$dist, "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("Log-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
  }

  invisible(x)
}

terms.crch <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.crch <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.crch <- function(object, model = c("location", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

fitted.crch <- function(object, type = c("location", "scale"), ...) object$fitted.values[[match.arg(type)]]

predict.crch <- function(object, newdata = NULL,
  type = c("response", "location", "scale", "quantile"), na.action = na.pass, at = 0.5, left = NULL, right = NULL, ...)
{
  type <- match.arg(type)
  ## response/location are synonymous
  if(type == "location") type <- "response"
  
  ## type quantile not available for user defined distributions
  if(type == "quantile" & identical(object$dist, "user defined"))
    stop("type quantile not available for user defined distributions")

  if(type == "quantile" & !is.null(newdata)){
      if(length(object$cens$left) > 1) {
        if(is.null(left)) stop("left has to be specified for non-constant left censoring")
        if(length(left) > 1 & length(left) != NROW(newdata)) stop("left must have length 1 or length of newdata")
      }
      if(length(object$cens$right) > 1) {
        if(is.null(right)) stop("right has to be specified for non-constant right censoring")
        if(length(right) > 1 & length(right) != NROW(newdata)) stop("right  must have length 1 or length of newdata")
      }
    }
  
  if(type == "quantile") {
    if(object$truncated) {
      qdist2 <- switch(object$dist, 
        "student"  = qtt, 
        "gaussian" = function(..., df) qtnorm(...), 
        "logistic" = function(..., df) qtlogis(...))
    } else {
      qdist2 <- switch(object$dist, 
        "student"  = qct, 
        "gaussian" = function(..., df) qcnorm(...), 
        "logistic" = function(..., df) qclogis(...))
    }
    
    if(is.null(left)) left <- object$cens$left
    if(is.null(right)) right <- object$cens$right

    qdist <- function(at, location, scale, df) {
      rval <- sapply(at, function(p) qdist2(p, location, scale, 
        df = df, left = left, right = right))
      if(length(at) > 1L) {
        if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
	        dimnames = list(unique(names(rval)), NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      } else {
        rval <- drop(rval)
      }
      rval 
    }
  }
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values$location
      },
      "scale" = {
        object$fitted.values$scale
      },
      "quantile" = {
        qdist(at, object$fitted.values$location, object$fitted.values$scale, object$df)
      },
    )
    return(rval)

  } else {
    tnam <- switch(type,
      "response" = "location",
      "scale" = "scale",
      "quantile" = "full"
    )

    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    newdata <- newdata[rownames(mf), , drop = FALSE]
    offset <- list(location = rep.int(0, nrow(mf)), scale = rep.int(0, nrow(mf)))

    if(type %in% c("response", "quantile")) {
      X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
      if(!is.null(object$call$offset)) offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
      if(!is.null(off.num <- attr(object$terms$location, "offset"))) {
        for(j in off.num) offset[[1L]] <- offset[[1L]] + eval(attr(object$terms$location, "variables")[[j + 1L]], newdata)
      }
    }
    if(type %in% c("scale", "quantile")) {
      Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)
      if(!is.null(off.num <- attr(object$terms$scale, "offset"))) {
        for(j in off.num) offset[[2L]] <- offset[[2L]] + eval(attr(object$terms$scale, "variables")[[j + 1L]], newdata)
      }
    }
   
    rval <- switch(type,
      "response" = {
        drop(X %*% object$coefficients$location + offset[[1L]])
      },
      "scale" = {
        object$link$scale$linkinv(drop(Z %*% object$coefficients$scale + offset[[2L]]))
      },
      "quantile" = {
        mu <- drop(X %*% object$coefficients$location + offset[[1L]])
        sigma <- object$link$scale$linkinv(drop(Z %*% object$coefficients$scale + offset[[2L]]))
        df <- object$df
        qdist(at, mu, sigma, df)       
      }
    )
    return(rval)

  }
}


coef.crch <- function(object, model = c("full", "location", "scale", "df"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = {
      cf$location
    },
    "scale" = {
      cf$scale
    },
    "df" = {
      cf$df
    },
    "full" = {
      nam <- c(names(cf$location), paste("(scale)", names(cf$scale), sep = "_"))
      if(length(cf$df)) nam <- c(nam, "(Log(df))")
      cf <- c(cf$location, cf$scale, cf$df)
      names(cf) <- nam
      cf
    }
  )
}


vcov.crch <- function(object, model = c("full", "location", "scale", "df"), ...) {
  vc <- object$vcov
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  l <- length(object$coefficients$df)

  model <-  match.arg(model)

  switch(model,
    "location" = {
      vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
    },
    "scale" = {
      vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$scale)
      vc
    },
    "df" = {
      vc <- vc[seq.int(length.out = l) + k + m, seq.int(length.out = l) + k + m, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$df)
      vc
    },
    "full" = {
      vc
    }
  )
}

logLik.crch <- function(object, ...) structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")

residuals.crch <- function(object, type = c("standardized", "pearson", "response", "quantile"), ...) {
  if(match.arg(type) == "response") {
    object$residuals 
  } else if (match.arg(type) == "quantile") {
    if(object$truncated) {
      pdist <- switch(object$dist, 
        "student"  = function(q, mean, sd) ptt(q, mean, sd, df = object$df), 
        "gaussian" = ptnorm, 
        "logistic" = ptlogis)
    } else {
      pdist <- switch(object$dist, 
        "student"  = function(q, mean, sd) pct(q, mean, sd, df = object$df), 
        "gaussian" = pcnorm, 
        "logistic" = pclogis)
    }
    qprob <- with(object, ifelse(residuals + fitted.values$location == cens$left,
      runif(n, 0, pdist(cens$left, fitted.values$location, fitted.values$scale)),
      ifelse(residuals + fitted.values$location == cens$right,
        runif(n, pdist(cens$right, fitted.values$location, fitted.values$scale), 1),
        pdist(residuals, 0, fitted.values$scale))))
    qnorm(qprob)
  } else object$residuals/object$fitted.values$scale
}



getSummary.crch <- function (obj, alpha = 0.05, ...) 
{
  cf <- summary(obj)$coefficients
  cval <- qnorm(1 - alpha/2)
  for (i in seq_along(cf)) cf[[i]] <- cbind(cf[[i]], 
      cf[[i]][, 1] - cval * cf[[i]][, 2],
      cf[[i]][, 1] + cval * cf[[i]][, 2])
  nam <- unique(unlist(lapply(cf, rownames)))
  acf <- array(dim = c(length(nam), 6, length(cf)), 
    dimnames = list(nam, c("est", "se", "stat", "p", "lwr", "upr"), names(cf)))
  for (i in seq_along(cf)) acf[rownames(cf[[i]]), , i] <- cf[[i]]

  return(list(
    coef = acf, 
    sumstat = c(
      N = obj$n, 
      logLik = as.vector(logLik(obj)), 
      AIC = AIC(obj), 
      BIC = BIC(obj)
    ), 
    contrasts = obj$contrasts, 
    xlevels = obj$xlevels, 
    call = obj$call))
}


update.crch <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}


