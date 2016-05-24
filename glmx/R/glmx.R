glmx <- function(formula, data, subset, na.action, weights, offset,
  family = negative.binomial, xlink = "log", control = glmx.control(...),
  model = TRUE, y = TRUE, x = FALSE, ...)
{  
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  Y <- model.response(mf, "any")
  X <- model.matrix(mt, mf)

  ## process response
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  n <- NROW(Y)
  if(n < 1) stop("empty model")

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  offset <- model.offset(mf)

  ## call the actual workhorse: glmx.fit()
  rval <- glmx.fit(X, Y, weights, offset, family, xlink, control)

  ## further model information
  rval$call <- cl
  rval$terms <- mt
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- attr(X, "contrasts")
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X

  class(rval) <- "glmx"  
  return(rval)
}

glmx.control <- function(profile = TRUE, nuisance = FALSE,
  start = NULL, xstart = NULL, hessian = TRUE,
  method = "BFGS", epsilon = 1e-8, maxit = c(500, 25), trace = FALSE,
  reltol = .Machine$double.eps^(1/1.2), ...)
{
  if(identical(hessian, TRUE)) hessian <- "optim"
  if(identical(hessian, FALSE)) hessian <- "none"
  hessian <- match.arg(hessian, c("none", "numDeriv", "optim"))

  maxit <- rep(maxit, length.out = 2L)
  trace <- rep(trace, length.out = 2L)
  glmcontrol <- glm.control(epsilon = epsilon, maxit = maxit[2L], trace = trace[2L])

  rval <- list(profile = profile, nuisance = nuisance,
    start = start, xstart = xstart, hessian = hessian, method = method,
    maxit = maxit[1L], trace = trace[1L], reltol = reltol,
    glm.control = glmcontrol)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- 1
  rval
}

glmx.fit <- function(x, y, weights = NULL, offset = NULL,
  family = negative.binomial, xlink = "log", control = glmx.control())
{
  ## process control options
  xctrl <- control
  profile <- control$profile
  nuisance <- control$nuisance
  glmstart <- control$start
  xstart <- control$xstart
  glmctrl <- control$glm.control
  method <- control$method
  hessian <- control$hessian
  xctrl$profile <- xctrl$nuisance <- xctrl$start <- xctrl$xstart <- xctrl$hessian <- xctrl$glm.control <- xctrl$method <- NULL
  
  ## process link for extra parameters
  if(!inherits(xlink, "link-glm")) xlink <- make.link(xlink)  

  ## starting values for extra parameters
  if(is.null(xstart)) xstart <- 0
  if(is.null(names(xstart))) {
    names(xstart) <- if(length(xstart) > 1L) {
      paste(names(formals(family))[1L], 1:length(xstart), sep = "")
    } else {
      names(formals(family))[1L]
    }
  }

  ## starting family
  family_start <- family(xlink$linkinv(xstart))

  ## response
  nobs <- n <- NROW(x)
  if(is.null(weights)) weights <- rep.int(1L, n)
  if(is.null(offset)) offset <- rep.int(0L, n)
  start <- etastart <- mustart <- NULL
  eval(family_start$initialize)

  ## update starting values for coefficients
  suppressWarnings(glmstart <- glm.fit(x, y, weights = weights, offset = offset,
    start = glmstart, control = glmctrl, family = family_start)$coefficients)

  ## regressors and parameters
  nobs <- sum(weights > 0L)
  k <- NCOL(x)
  q <- length(xstart)

  if(is.null(family_start$dispersion)) family_start$dispersion <- family_start$family %in% c("gaussian", "Gamma", "inverse.gaussian")
  if(family_start$dispersion) {
    dispersion <- function(wresiduals, wweights) sum(wresiduals^2, na.rm = TRUE)/sum(wweights, na.rm = TRUE)
    dpar <- 1
  } else {
    dispersion <- function(wresiduals, wweights) 1
    dpar <- 0
  }
  df <- k + q + dpar

  ## objective function
  profile_loglik <- function(par) {
    suppressWarnings(gm <- glm.fit(x, y, weights = weights, offset = offset,
      start = glmstart, control = glmctrl, family = family(xlink$linkinv(par))))
    aic <- gm$aic
    aic/2 - dpar - length(glmstart)
  }

  profile_grad <- if(is.null(family_start$loglik.extra)) NULL else function(par) {
    gamma <- par
    extra <- xlink$linkinv(gamma)
    f <- family(extra)
    suppressWarnings(beta <- glm.fit(x, y, weights = weights, offset = offset,
      start = glmstart, control = glmctrl, family = f)$coefficients)
    eta <- drop(x %*% beta + offset)
    fog <- f$mu.eta(eta)
    mu <- f$linkinv(eta)    
    varmu <- f$variance(mu)    
    phi <- dispersion((y - mu) / varmu, fog)
    ggamma <- matrix(f$loglik.extra(y, mu, extra) * xlink$mu.eta(gamma), nrow = length(mu))
    -colSums(ggamma/phi)
  }

  full_loglik <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    f <- family(xlink$linkinv(gamma))
    mu <- f$linkinv(drop(x %*% beta + offset))    
    dev <- sum(f$dev.resids(y, mu, weights))
    f$aic(y, n, mu, weights, dev)/2 - dpar
  }

  full_grad <- if(is.null(family_start$loglik.extra)) NULL else function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    extra <- xlink$linkinv(gamma)
    f <- family(extra)
    eta <- drop(x %*% beta + offset)
    fog <- f$mu.eta(eta)
    mu <- f$linkinv(eta)    
    varmu <- f$variance(mu)    
    phi <- dispersion((y - mu) / varmu, fog)
    gbeta <- sqrt(weights) * ((y - mu) / varmu) * fog
    ggamma <- matrix(f$loglik.extra(y, mu, extra) * xlink$mu.eta(gamma), nrow = length(mu))
    -colSums(cbind(gbeta * x, ggamma)/phi)
  }

  if(profile) {
    ## optimize profile likelihood first
    opt1 <- optim(par = xstart, fn = profile_loglik, gr = profile_grad,
      method = method, hessian = FALSE, control = xctrl)
    if(opt1$convergence > 0L) warning("optimization of profile likelihood failed to converge")

    ## extract optimal parameters
    xpar <- opt1$par

    ## optimal coefficients at chosen extra parameters
    cf <- glm.fit(x, y, weights = weights, offset = offset,
      start = glmstart, control = glmctrl, family = family(xlink$linkinv(xpar)))$coefficients
  } else {
    ## use starting values of extra parameters
    opt1 <- NULL
    xpar <- xstart
    cf <- glm.fit(x, y, weights = weights, offset = offset,
      start = glmstart, control = glmctrl, family = family(xlink$linkinv(xpar)))$coefficients
  }

  ## optimize full likelihood
  opt2 <- optim(par = c(cf, xpar), fn = full_loglik, gr = full_grad,
    method = method, hessian = hessian == "optim", control = xctrl)
  if(opt2$convergence > 0L) warning("optimization of full likelihood failed to converge")

  ## extract fitted values/parameters
  beta <- as.vector(opt2$par[1:k])
  gamma <- as.vector(opt2$par[-(1:k)])
  f <- family(xlink$linkinv(gamma))
  eta <- drop(x %*% beta + offset)
  mu <- f$linkinv(eta)
  dev <- f$dev.resids(y, mu, weights)
  phi <- dispersion((y - mu) / f$variance(mu), f$mu.eta(eta))

  ## names
  if(!is.null(colnames(x))) {
    names(beta) <- colnames(x)
    names(gamma) <- names(xstart)
    if(xlink$name != "identity") names(gamma) <- paste(xlink$name, "(", names(gamma), ")", sep = "")
  }

  if(hessian != "none") {
    vc <- if(hessian == "numDeriv") {
      numDeriv::hessian(full_loglik, c(beta, gamma))
    } else {
      as.matrix(opt2$hessian)
    }
    vc <- solve(vc)
    rownames(vc) <- colnames(vc) <- c(names(beta), names(gamma))
  } else {
    vc <- NULL
  }

  ## set up return value
  n <- NROW(x)
  rval <- list(  
    coefficients = list(glm = beta, extra = gamma),
    residuals = dev,
    fitted.values = mu,
    optim = list(profile = opt1, full = opt2),
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, list(mean = rep(0, n), scale = rep(0, n)))) NULL else offset,
    n = n,
    nobs = nobs,
    df = df,
    loglik = -opt2$value,
    dispersion = phi,
    vcov = vc,
    family = list(glm = f, extra = family),
    xlink = xlink,
    control = control,
    converged = opt2$convergence < 1
  )
  class(rval) <- "glmx"

  return(rval)
}

formula.glmx <- function(x, ...) x$call$formula

coef.glmx <- function(object, model = NULL, ...) {
  if(is.null(model)) {
    model <- if(object$control$nuisance) "glm" else "full"
  }
  if(model == "x") model <- "extra"
  switch(model,
    "full" = c(object$coefficients$glm, object$coefficients$extra),
    "glm" = object$coefficients$glm,
    "extra" = object$coefficients$extra
  )
}

vcov.glmx <- function(object, model = NULL, ...) {
  if(is.null(model)) {
    model <- if(object$control$nuisance) "glm" else "full"
  }
  if(model == "x") model <- "extra"
  k <- length(object$coefficients$glm)
  switch(model,
    "full" = object$vcov,
    "glm" = object$vcov[1L:k, 1L:k, drop = FALSE],
    "extra" = object$vcov[-(1L:k), -(1L:k), drop = FALSE]
  )
}

nobs.glmx <- function(object, ...) object$nobs

logLik.glmx <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

print.glmx <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(coef(x))) {
      cat(paste("Coefficients (", x$family$glm$family, " model with ", x$family$glm$link, " link):\n", sep = ""))
      print.default(format(x$coefficients$glm, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
      if(!x$control$nuisance) {
        cat(paste("Extra parameters (with ", x$xlink$name, " link):\n", sep = ""))
        print.default(format(x$coefficients$extra, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      }
    }
    else cat("No coefficients\n\n")
  }
  
  invisible(x)
}

summary.glmx <- function(object, vcov. = NULL, type = "deviance", ...)
{
  ## residuals
  type <- match.arg(type, "deviance") #FIXME# c("deviance", "pearson", "response"))
  #FIXME# object$residuals <- residuals(object, type = type)
  object$residuals.type <- type
  
  ## extend coefficient table
  cf <- coef(object)
  se <- if(is.null(vcov.)) sqrt(diag(object$vcov)) else sqrt(diag(vcov.(object)))
  if(!is.null(names(se))) se <- se[names(cf)]
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  k <- length(object$coefficients$glm)
  object$coefficients <- list(glm = cf[1:k, , drop = FALSE], extra = cf[-(1:k), , drop = FALSE])
  
  ## number of iterations
  object$iterations <- c(
    "profile" = if(is.null(object$optim$profile)) 0 else object$optim$profile$counts["gradient"],
    "full" = object$optim$full$counts["gradient"]
  )
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.glmx"
  object
}

print.summary.glmx <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    types <- c("deviance", "pearson", "response")
    Types <- c("Deviance residuals", "Pearson residuals", "Raw response residuals")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  
    cat(paste("\nCoefficients (", x$family$glm$family, " model with ", x$family$glm$link, " link):\n", sep = ""))
    printCoefmat(x$coefficients$glm, digits = digits, signif.legend = FALSE)
  
    if(!x$control$nuisance) {
      cat(paste("\nExtra parameters (with ", x$xlink$name, " link):\n", sep = ""))
        printCoefmat(x$coefficients$extra, digits = digits, signif.legend = FALSE)
    }
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[,4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    cat("\nDispersion:", if(x$dispersion == 1) "1" else formatC(x$dispersion, digits = digits))
    cat(paste("\nNumber of iterations in", x$control$method, "optimization:",
      x$iterations[1L], "(profile)", x$iterations[2L], "(full)\n"))
  }
  
  invisible(x)
}

predict.glmx <- function(object, newdata = NULL,
  type = c("response", "link"), na.action = na.pass, ...) 
{
  type <- match.arg(type)
  f <- object$family$glm
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = object$fitted.values,
      "link" = f$linkfun(object$fitted.values)
    )
    names(rval) <- names(object$fitted.values)
    return(rval)

  } else {

    ## model frame and matrix
    mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)

    ## offset
    newdata <- newdata[rownames(mf), , drop = FALSE]
    offset <- rep.int(0, nrow(mf))
    if(!is.null(object$call$offset)) offset <- offset + eval(object$call$offset, newdata)
    if(!is.null(off.num <- attr(object$terms, "offset"))) {
      for(j in off.num) offset <- offset + eval(attr(object$terms, "variables")[[j + 1L]], newdata)
    }

    ## linear predictor
    eta <- drop(X %*% object$coefficients$glm + offset)

    rval <- switch(type,    
      "response" = f$linkinv(eta),      
      "link" = eta
    )
    names(rval) <- rownames(mf)
    return(rval)

  }
}
