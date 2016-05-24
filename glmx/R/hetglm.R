hetglm <- function(formula, data, subset, na.action, weights, offset,
                   family = binomial(link = "probit"),
                   link.scale = c("log", "sqrt", "identity"),
 		   control = hetglm.control(...),
		   model = TRUE, y = TRUE, x = FALSE, ...)
{  
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  keep_y <- y
  keep_x <- x
  
  ## scale link
  if(is.character(link.scale)) link.scale <- match.arg(link.scale)

  ## family
  if(is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2] < 2L) {
    formula <- as.Formula(formula(formula), formula(formula, lhs = 0L))
  } else {
    if(length(formula)[2] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  attr(mtZ, "intercept") <- 1L
  Y <- model.response(mf, "any")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)[, -1L, drop = FALSE]

  ## process response
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  n <- NROW(Y)
  if(n < 1) stop("empty model")

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1L) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## in mean part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## in scale part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for mean)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(mean = offsetX, scale = offsetZ)

  ## initialize family (essentially: process response)
  nobs <- n
  y <- Y
  start <- etastart <- mustart <- NULL
  eval(family$initialize)
  Y <- y  

  ## call the actual workhorse: hetglm.fit()
  rval <- hetglm.fit(X, Y, Z, weights, offset, family, link.scale, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, scale = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(keep_y) rval$y <- Y
  if(keep_x) rval$x <- list(mean = X, scale = Z)

  class(rval) <- "hetglm"  
  return(rval)
}

hetglm.control <- function(method = "nlminb", maxit = 1000, hessian = FALSE, trace = FALSE, start = NULL, ...)
{
  if(method == "nlminb") {
    rval <- list(method = method, trace = as.numeric(trace), start = start)
    rval <- c(rval, list(...))
    if(is.null(rval$iter.max)) rval$iter.max <- maxit  
    if(is.null(rval$eval.max)) rval$eval.max <- max(200, rval$iter.max)
  } else {
    rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start)
    rval <- c(rval, list(...))
    if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
    rval$fnscale <- 1
  }
  rval
}

hetglm.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
  family = binomial(), link.scale = "log", control = hetglm.control())
{
  ## response
  nobs <- n <- NROW(x)
  if(is.null(weights)) weights <- rep.int(1L, n)
  if(is.null(offset)) offset <- list(mean = rep.int(0, n), scale = rep.int(0, n))
  start <- etastart <- mustart <- NULL
  eval(family$initialize)

  ## regressors
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(z)) z <- x[, -1L, drop = FALSE]
  m <- NCOL(z)

  ## account for estimated dispersion
  if(is.null(family$dispersion)) family$dispersion <- family$family %in% c("gaussian", "Gamma", "inverse.gaussian")
  if(family$dispersion) {
    dispersion <- function(wresiduals, wweights) sum(wresiduals^2, na.rm = TRUE)/sum(wweights, na.rm = TRUE)
    dpar <- 1
  } else {
    dispersion <- function(wresiduals, wweights) 1
    dpar <- 0
  }
  ## FIXME: Not clear whether dispersion is handled consistently in
  ## "glm". In summary.glm() the dispersion is only fixed for "binomial"
  ## and "poisson". In logLik.glm() an additional parameter is accounted for
  ## only for "gaussian", "Gamma", "inverse.gaussian".
  ## Maybe suggest adding disperion = FALSE/TRUE to family objects?

  ## link processing
  linkinv <- family$linkinv
  linkobj <- family[c("linkfun", "linkinv", "mu.eta", "valideta", "link")]
  names(linkobj)[5] <- "name"
  class(linkobj) <- "link-glm"
  scale_linkstr <- link.scale
  variance <- family$variance
  mu.eta <- family$mu.eta
  dev.resids <- family$dev.resids
  aic <- family$aic
  scale_linkobj <- make.link(scale_linkstr)
  scale_linkfun <- scale_linkobj$linkfun
  scale_linkinv <- scale_linkobj$linkinv
  scale_mu.eta <- scale_linkobj$mu.eta

  ## control parameters
  method <- control$method
  start <- control$start
  if(method == "nlminb") {
    hessian <- FALSE
  } else {
    hessian <- control$hessian
    control$hessian <- NULL
  }
  control <- control[-which(names(control) %in% c("method", "start"))]

  ## null model and default starting values
  nullreg <- glm.fit(x = x, y = y, weights = weights, offset = offset$mean, family = family)
  if(is.null(start)) {
    start <- list(
      mean = nullreg$coefficients,
      scale = rep.int(0, m)
    )
  }
  if(is.list(start)) start <- do.call("c", start)

  ## objective function and gradient
  loglikfun <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    scale_eta <- scale_linkfun(1) + drop(z %*% gamma + offset$scale)
    scale <- scale_linkinv(scale_eta)
    eta <- drop(x %*% beta + offset$mean) / scale
    mu <- linkinv(eta)
    
    dev <- sum(dev.resids(y, mu, weights))
    aic(y, n, mu, weights, dev)/2 - dpar
  }
  gradfun <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    scale_eta <- scale_linkfun(1) + drop(z %*% gamma + offset$scale)
    scale <- scale_linkinv(scale_eta)
    eta <- drop(x %*% beta + offset$mean) / scale
    fog <- mu.eta(eta) / scale ## aka working weights
    mu <- linkinv(eta)
    varmu <- variance(mu)
    phi <- dispersion((y - mu) / varmu, fog)
 
    gbeta <- sqrt(weights) * ((y - mu) / varmu) * fog
    ggamma <- - gbeta * eta * scale_mu.eta(scale_eta)
    -colSums(cbind(gbeta * x, ggamma * z)/phi)
  }
  hessfun <- function(par, inverse = FALSE) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    scale_eta <- scale_linkfun(1) + drop(z %*% gamma + offset$scale)
    scale <- scale_linkinv(scale_eta)
    eta <- drop(x %*% beta + offset$mean) / scale
    fog <- mu.eta(eta) / scale
    mu <- linkinv(eta)
    varmu <- variance(mu)   
    phi <- dispersion((y - mu) / varmu, fog)
 
    Hbeta <- sqrt(weights) * (1 / sqrt(varmu)) * fog
    Hgamma <- - Hbeta * eta * scale_mu.eta(scale_eta)
    if(!inverse) {
      crossprod(cbind(Hbeta * x, Hgamma * z))/phi
    } else {  ## better than: solve(crossprod(...))
      chol2inv(qr.R(qr(cbind(Hbeta * x, Hgamma * z)))) * phi
    }
  }

  ## optimize likelihood  
  if(method == "nlminb") {
    opt <- nlminb(start = start, objective = loglikfun, gradient = gradfun,
      hessian = hessfun, control = control)  
    if(opt$convergence > 0) warning(paste("optimization failed to converge with message:",
      opt$message))
  } else {
    opt <- optim(par = start, fn = loglikfun, gr = gradfun,
      method = method, hessian = hessian, control = control)
    if(opt$convergence > 0) warning("optimization failed to converge")
  
  }


  ## extract fitted values/parameters
  vc <- if(hessian) solve(as.matrix(opt$hessian)) else hessfun(opt$par, inverse = TRUE)
  beta <- as.vector(opt$par[1:k])
  gamma <- as.vector(opt$par[-(1:k)])
  eta <- drop(x %*% beta + offset$mean)
  scale_eta <- scale_linkfun(1) + drop(z %*% gamma + offset$scale)
  scale <- scale_linkinv(scale_eta)
  mu <- linkinv(eta / scale)
  phi <- dispersion((y - mu) / variance(mu), mu.eta(eta) / scale)
  nobs <- sum(weights > 0L)  


  ## names
  if(!is.null(colnames(x))) {
    names(beta) <- colnames(x)
    if(!is.null(colnames(z))) names(gamma) <- colnames(z)
    rownames(vc) <- colnames(vc) <- c(colnames(x),
      if(m > 0L) paste("(scale)", colnames(z), sep = "_") else NULL)
  }

  ## set up return value
  rval <- list(  
    coefficients = list(mean = beta, scale = gamma),
    residuals = y - mu,
    fitted.values = structure(mu, .Names = names(y)),
    optim = opt,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, list(mean = rep(0, n), scale = rep(0, n)))) NULL else offset,
    n = n,
    nobs = nobs,
    df.null = nobs - k,
    df.residual = nobs - k - m,
    loglik = -loglikfun(opt$par),
    loglik.null = -loglikfun(c(nullreg$coefficients, rep.int(0, m))),
    dispersion = phi,
    vcov = vc,
    family = family,
    link = list(mean = linkobj, scale = scale_linkobj),
    converged = opt$convergence < 1    
  )

  return(rval)
}

print.hetglm <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(coef(x))) {
      cat(paste("Coefficients (", x$family$family, " model with ", x$link$mean$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
      cat(paste("Latent scale model coefficients (with ", x$link$scale$name, " link):\n", sep = ""))
      if(length(x$coefficients$scale)) {
        print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      } else {
        cat("None (constant scale = 1).")
      }
      cat("\n")
    }
    else cat("No coefficients\n\n")
  }
  
  invisible(x)
}

summary.hetglm <- function(object, vcov. = NULL, type = "deviance", ...)
{
  ## residuals
  type <- match.arg(type, c("deviance", "pearson", "response"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type
  
  ## extend coefficient table
  k <- length(object$coefficients$mean)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- if(is.null(vcov.)) sqrt(diag(object$vcov)) else sqrt(diag(vcov.(object)))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[1:k, , drop = FALSE], scale = cf[-(1:k), , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf
  
  ## LR test against null model
  object$lrtest <- if(length(object$coefficients$scale)) {
    c("LR" = -2 * (object$loglik.null - object$loglik),
      "Df" = object$df.null - object$df.resid,
      "p-value" = pchisq(-2 * (object$loglik.null - object$loglik),
        object$df.null - object$df.resid, lower.tail = FALSE))
  } else {
    rep(NA, 3)
  }
  
  ## number of iterations
  iter <- if(object$method == "nlminb") "evaluations" else "counts"
  object$iterations <- as.vector(na.omit(object$optim[[iter]]["gradient"]))
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.hetglm"
  object
}

print.summary.hetglm <- function(x, digits = max(3, getOption("digits") - 3), ...)
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
  
    cat(paste("\nCoefficients (", x$family$family, " model with ", x$link$mean$name, " link):\n", sep = ""))
    printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
  
    cat(paste("\nLatent scale model coefficients (with ", x$link$scale$name, " link):\n", sep = ""))
    if(length(x$coefficients$scale)) {
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else {
      cat("None (constant scale = 1).\n")
    }
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[,4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$lrtest[1])) cat("\nLR test for homoskedasticity:",
      formatC(x$lrtest[1], digits = digits), "on", x$lrtest[2], "Df, p-value:", format.pval(x$lrtest[3], digits = digits))
    cat("\nDispersion:", if(x$dispersion == 1L) "1" else formatC(x$dispersion, digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
  }
  
  invisible(x)
}

predict.hetglm <- function(object, newdata = NULL,
  type = c("response", "link", "scale"), na.action = na.pass, ...) 
{
  type <- match.arg(type)
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values
      },
      "link" = {
        object$link$mean$linkfun(object$fitted.values)
      },      
      "scale" = {
        gamma <- object$coefficients$scale
        z <- if(is.null(object$x)) model.matrix(object, model = "scale") else object$x$scale
	offsetz <- if(is.null(object$offset)) 0 else offset$scale
	object$link$scale$linkinv(object$link$scale$linkfun(1) + drop(z %*% gamma + offsetz))
      }
    )
    names(rval) <- names(object$fitted.values)
    return(rval)

  } else {

    ## model frame and matrices
    mf <- model.frame(delete.response(object$terms[["full"]]), newdata, na.action = na.action, xlev = object$levels[["full"]])
    X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
    Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)[, -1L, drop = FALSE]

    ## offsets
    newdata <- newdata[rownames(mf), , drop = FALSE]
    offset <- list(mean = rep.int(0, nrow(mf)), scale = rep.int(0, nrow(mf)))
    if(!is.null(object$call$offset)) offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
    if(!is.null(off.num <- attr(object$terms$mean, "offset"))) {
      for(j in off.num) offset[[1L]] <- offset[[1L]] + eval(attr(object$terms$mean, "variables")[[j + 1L]], newdata)
    }
    if(!is.null(off.num <- attr(object$terms$scale, "offset"))) {
      for(j in off.num) offset[[2L]] <- offset[[2L]] + eval(attr(object$terms$scale, "variables")[[j + 1L]], newdata)
    }

    ## linear predictors
    eta_mean <- drop(X %*% object$coefficients$mean + offset$mean)
    pred_scale <- object$link$scale$linkinv(object$link$scale$linkfun(1) + drop(Z %*% object$coefficients$scale + offset$scale))

    rval <- switch(type,    
      "response" = {
        object$link$mean$linkinv(eta_mean / pred_scale)
      },      
      "link" = {
        eta_mean / pred_scale
      },      
      "scale" = {
        pred_scale
      }
    )
    names(rval) <- rownames(mf)
    return(rval)

  }
}

coef.hetglm <- function(object, model = c("full", "mean", "scale"), ...) {
  cf <- object$coefficients

  model <-  match.arg(model)  
  switch(model,
    "mean" = {
      cf$mean
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      cf <- c(cf$mean, cf$scale)
      names(cf) <- colnames(object$vcov)
      cf
    }
  )
}

vcov.hetglm <- function(object, model = c("full", "mean", "scale"), ...) {
  vc <- object$vcov
  k <- length(object$coefficients$mean)

  model <- match.arg(model)
  switch(model,
    "mean" = {
      vc[1:k, 1:k, drop = FALSE]
    },
    "scale" = {
      vc <- vc[-(1:k), -(1:k), drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$scale)
      vc
    },
    "full" = {
      vc
    }
  )
}

bread.hetglm <- function(x, ...) {
  vcov(x) * x$nobs
}

estfun.hetglm <- function(x, ...)
{
  ## extract response y and regressors X and Z
  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  xmat <- if(is.null(x$x)) model.matrix(x, model = "mean") else x$x$mean
  zmat <- if(is.null(x$x)) model.matrix(x, model = "scale") else x$x$scale
  offset <- if(is.null(x$offset)) list(mean = 0, scale = 0) else x$offset
  wts <- weights(x)
  if(is.null(wts)) wts <- 1

  if(x$family$dispersion) {
    dispersion <- function(wresiduals, wweights) sum(wresiduals^2, na.rm = TRUE)/sum(wweights, na.rm = TRUE)
  } else {
    dispersion <- function(wresiduals, wweights) 1
  }

  beta <- x$coefficients$mean
  gamma <- x$coefficients$scale
  scale_eta <- x$link$scale$linkfun(1) + drop(zmat %*% gamma + offset$scale)
  scale <- x$link$scale$linkinv(scale_eta)
  eta <- drop(xmat %*% beta + offset$mean) / scale
  mu <- x$link$mean$linkinv(eta)
  wt <- x$link$mean$mu.eta(eta) / scale
  resid <- (y - mu) / x$family$variance(mu)
  phi <- dispersion(resid, wt)

  gbeta <- resid * x$link$mean$mu.eta(eta) / scale
  ggamma <- - gbeta * eta * x$link$scale$mu.eta(scale_eta)

  rval <- cbind(wts * gbeta * xmat/phi, wts * ggamma * zmat/phi)
  rownames(rval) <- names(y)
  colnames(rval) <- colnames(x$vcov)

  attr(rval, "assign") <- NULL
  return(rval)
}

coeftest.hetglm <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)  

logLik.hetglm <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)) + object$family$dispersion, class = "logLik")
}

terms.hetglm <- function(x, model = c("mean", "scale"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.hetglm <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.Formula(formula$formula)) formula$call$formula <- formula$formula <-
    formula(formula$formula, collapse = TRUE)
  formula$terms <- formula$terms$full
  NextMethod()
}

model.matrix.hetglm <- function(object, model = c("mean", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) {
    object$x[[model]]
  } else {
    mm <- model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
    if(model == "scale") mm[, -1L, drop = FALSE] else mm
  }
  return(rval)
}

residuals.hetglm <- function(object,
  type = c("deviance", "pearson", "response"), ...)
{
  ## desired type
  type <- match.arg(type)

  ## extract fitted information
  res <- object$residuals
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  mu <- fitted(object)
  wts <- weights(object)
  if(is.null(wts)) wts <- rep.int(1, length(mu))
  
  res <- switch(type,  
    "deviance" = {
      sign(res) * sqrt(object$family$dev.resids(y, mu, wts))
    },
    "pearson" = {
      sqrt(wts) * res / sqrt(object$family$variance(mu))
    },
    "response" = {
      sqrt(wts) * res
    }
  )
  names(res) <- names(mu)

  return(res)
}

update.hetglm <- function (object, formula., ..., evaluate = TRUE)
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
