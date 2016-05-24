heavyLm <-
function(formula, data, family = Student(df = 4), subset, na.action,
    control = heavy.control(), model = TRUE, x = FALSE, y = FALSE,
    contrasts = NULL)
{
    ret.x <- x
    ret.y <- y
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    is.multivariate <- is.matrix(y)
    x <- model.matrix(Terms, mf, contrasts)
    xnames <- dimnames(x)[[2]]
    if (is.multivariate)
      ynames <- dimnames(y)[[2]]
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    
    ## initial estimates
    fit <- lsfit(x, y, intercept = FALSE)[1:2]
    res <- fit$residuals
    cf <- fit$coefficients
    if (is.multivariate)
      Sigma <- crossprod(res) / n
    else
      sigma2 <- sum(res^2) / n
    
    ## extract family info
    if (!inherits(family, "heavy.family"))
      stop("Use only with 'heavy.family' objects")
    if (is.null(family$family))
      stop("'family' not recognized")
    kind <- family$which
    if ((kind < 0) || (kind > 4))
      stop("not valid 'family' object")
    settings <- c(kind, family$npars, unlist(family$pars))
    
    ## set control values
    if (missing(control))
      control <- heavy.control()
    if (!control$algorithm)
      control$ncycles <- 1
    ctrl <- unlist(control)[1:4]
    ctrl <- c(ctrl, 0)
    
    ## Call fitter
    now <- proc.time()
    if (is.multivariate) {
      ny <- ncol(y)
      dx <- c(dx, ny)
      #storage.mode(y) <- "double"
      #storage.mode(x) <- "double"
      #storage.mode(res) <- "double"
      fit <- .C("mlm_fit",
                y = as.double(y),
                x = as.double(x),
                dims = as.integer(dx),
                settings = as.double(settings),
                coefficients = as.double(cf),
                Sigma = as.double(Sigma),
                fitted = double(n * ny),
                residuals = as.double(res),
                distances = double(n),
                weights = as.double(rep(1, n)),
                logLik = double(1),
                acov = double(p^2),
                control = as.double(ctrl))
    }
    else {
      fit <- .C("lm_fit",
                y = as.double(y),
                x = as.double(x),
                dims = as.integer(dx),
                settings = as.double(settings),
                coefficients = as.double(cf),
                sigma2 = as.double(sigma2),
                fitted = double(n),
                residuals = as.double(res),
                distances = double(n),
                weights = as.double(rep(1, n)),
                logLik = double(1),
                acov = double(p^2),
                control = as.double(ctrl))
    }
    speed <- proc.time() - now
    
    ## creating the output object
    out <- list(call = Call,
                dims = dx,
                family = family,
                settings = fit$settings)
    if (is.multivariate) {
      out$coefficients <- matrix(fit$coefficients, ncol = ny)
      out$Sigma <- matrix(fit$Sigma, ncol = ny)
      out$fitted.values <- matrix(fit$fitted, ncol = ny)
      out$residuals <- matrix(fit$residuals, ncol = ny)
    } else {
      out$coefficients <- fit$coefficients
      out$sigma2 <- fit$sigma2
      out$fitted.values <- fit$fitted
      out$residuals <- fit$residuals
    }
    out$logLik <- fit$logLik
    out$numIter <- fit$control[5]
    out$control <- control
    out$weights <- fit$weights
    out$distances <- fit$distances
    out$acov <- matrix(fit$acov, ncol = p)
    out$speed = speed
    out$converged = FALSE
    if (!control$fix.shape) {
      if ((kind > 1) && (kind < 4)) {
        df <- signif(out$settings[3], 6)
        out$family$call <- call(out$family$family, df = df)
      }
    }
    if (out$numIter < control$maxIter)
      out$converged <- TRUE
    if (is.multivariate) {
      dimnames(out$coefficients)[[1]] <- xnames
      dimnames(out$coefficients)[[2]] <- ynames
      dimnames(out$Sigma)[[1]] <- ynames
      dimnames(out$Sigma)[[2]] <- ynames
      out$acov <- kronecker(out$acov, out$Sigma)
    } else
      names(out$coefficients) <- xnames
    out$na.action <- attr(mf, "na.action")
    out$contrasts <- attr(x, "contrasts")
    out$xlevels <- .getXlevels(Terms, mf)
    out$terms <- Terms
    if (model)
      out$model <- mf
    if (ret.y)
      out$y <- y
    if (ret.x)
      out$x <- x
    if (is.multivariate)
      class(out) <- "heavyMLM"
    else
      class(out) <- "heavyLm"
    out
}

print.heavyLm <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded\n")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$sigma2), "\n")
  invisible(x)
}

print.heavyMLM <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded\n")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\n")
  invisible(x)
}

summary.heavyLm <-
function (object, ...)
{
  z <- object
  se <- sqrt(diag(z$acov))
  est <- z$coefficients
  zval <- est / se
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$family <- z$family
  ans$logLik <- z$logLik
  ans$sigma2 <- z$sigma2
  ans$residuals <- z$residuals
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients),
        c("Estimate", "Std.Error", "Z value", "p-value"))
  ans$correlation <- z$acov / outer(se, se)
  dimnames(ans$correlation) <- dimnames(ans$coefficients)[c(1,1)]
  class(ans) <- "summary.heavyLm"
  ans
}

summary.heavyMLM <-
function (object, ...)
{
  z <- object
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$family <- z$family
  ans$logLik <- z$logLik
  ans$coefficients <- z$coefficients
  ans$Sigma <- z$Sigma
  ans$residuals <- z$residuals
  ans$acov <- z$acov
  ans$correlation <- cov2cor(z$acov)
  class(ans) <- "summary.heavyMLM"
  ans
}

print.summary.heavyLm <-
function(x, digits = 4, ...)
{
  cat("Linear model under heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  resid <- x$residuals
  nobs <- x$dims[1]
  p <- x$dims[2]
  rdf <- nobs - p
  if (rdf > 5) {
    cat("\nResiduals:\n")
		rq <- quantile(resid)
		names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
		print(rq, digits = digits, ...)
	}
	else if(rdf > 0) {
	 cat("\nResiduals:\n")
	 print(resid, digits = digits, ...)
  }
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$sigma2))
  cat("\nLog-likelihood:", format(x$logLik), "on", p + 1, "degrees of freedom\n")
  invisible(x)
}

print.summary.heavyMLM <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Multivariate regression under heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  nobs <- x$dims[1]
  p <- x$dims[2]
  rdf <- nobs - p
  ny <- x$dims[3]
  npars <- ny * p + ny * (ny + 1) / 2
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  print.symmetric(x$Sigma, digits = digits)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nLog-likelihood:", format(x$logLik), "on", npars, "degrees of freedom\n")
  invisible(x)
}
