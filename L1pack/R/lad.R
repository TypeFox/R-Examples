lad <- 
function(formula, data, method = c("BR", "EM"), subset, na.action, 
  control = NULL, model = TRUE, x = FALSE, y = FALSE, contrasts = NULL)
{
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  xnames <- dimnames(x)[[2]]
  dx <- dim(x)
  n  <- dx[1]

  ## set control values 
  if (is.null(control))
    control <- l1pack.control()
  method <- match.arg(method)
  choice <- switch(method, "BR" = 0, "EM" = 1)
  ctrl  <- unlist(control)
  nctrl <- names(control)
  ctrl  <- c(ctrl, choice, 0)

  ## initial estimates
  fit <- lsfit(x, y, intercept = FALSE)[1:2]
  res <- fit$residuals
  cf <- fit$coefficients

  ## Call fitter
  now <- proc.time()
  fit <- .C("lad",
            y = as.double(y),
            x = as.double(x),
            dims = as.integer(dx),
            coefficients = as.double(cf),
            scale  = as.double(0),
            fitted = double(n),
            resid  = as.double(res),
            weights = as.double(rep(1, n)),
            control = as.double(ctrl),
            sad = as.double(0),
            logLik = as.double(0))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              dims = fit$dims,
              coefficients = fit$coefficients,
              scale = fit$scale,
              minimum = fit$sad,
              fitted.values = fit$fitted,
              residuals = fit$resid,
              numIter = fit$control[4],
              control = fit$control,
              weights = fit$weights,
              logLik = fit$logLik,
              speed = speed,
              converged = FALSE)
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  names(out$control) <- c(nctrl, "method", "numIter")
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
  class(out) <- "lad"
  out
}

print.lad <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$scale), "\n")
  invisible(x)
}
