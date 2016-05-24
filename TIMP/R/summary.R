summary.timp.optim <- function (object, currModel, ...)
{
  param <- object$par
  pnames <- names(param)
  ih <- solve(object$hessian)
  nnonlin <- length(param)
  p <- nnonlin + currModel@nclp
  ## have to get the deviance
  dev <- currModel@fit@rss
  ## and the number of obs
  obs<-0
  for(i in 1:length(currModel@data)) 
    obs <- obs + length(currModel@data[[i]]@psi.df) 
  rdf <- obs - nnonlin 
  resvar <- dev / rdf
  se <- sqrt(diag(ih) * resvar)
  names(se) <- pnames
  tval <- param/se
  param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans <- list(sigma = sqrt(dev/rdf),
              df = c(p, rdf), cov.unscaled = ih,
              coefficients = param, nnonlin = nnonlin,
              nclp = currModel@nclp)# could return more but don't bother now
  class(ans) <- "summary.timp.optim"
  ans
}
summary.timp.nls.lm <- function (object, currModel, ...)
{
  class(object) <- "nls.lm" # so that methods for 'nls.lm' work
  param <- coef(object)
  pnames <- names(param)
  ih <- solve(object$hessian)
  nnonlin <- length(param)
  p <- nnonlin + currModel@nclp
  rdf <- length(object$fvec) - nnonlin 
  resvar <- deviance(object) / rdf
  se <- sqrt(diag(ih) * resvar)
  names(se) <- pnames
  tval <- param/se
  
  param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans <- list(residuals = object$fvec, sigma = sqrt(object$deviance/rdf),
              df = c(p, rdf), cov.unscaled = ih,
              info = object$info, niter = object$niter,
              stopmess = object$message, 
              coefficients = param, nnonlin = nnonlin,
              nclp = currModel@nclp)
  class(ans) <- "summary.timp.nls.lm"
  ans
}
summary.timp.nls <-
    function (object, currModel, correlation = FALSE,
              symbolic.cor = FALSE, ...)
{
    class(object) <- "nls" # so that methods for 'nls' work
    r <- as.vector(object$m$resid()) # These are weighted residuals.
    w <- object$weights
    n <- if (!is.null(w)) sum(w > 0) else length(r)
    param <- coef(object)
    pnames <- names(param)
    nnonlin <- length(param)
    p <- nnonlin + currModel@nclp
    rdf <- n - nnonlin
    resvar <- if(rdf <= 0) NaN else deviance(object)/rdf
    XtXinv <- chol2inv(object$m$Rmat())
    dimnames(XtXinv) <- list(pnames, pnames)
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(formula = formula(object), residuals = r, sigma = sqrt(resvar),
                df = c(p, rdf), cov.unscaled = XtXinv,
                call = object$call,
                convInfo = object$convInfo,
                control = object$control,
                na.action = object$na.action,
                coefficients = param,
                parameters = param, nnonlin = nnonlin,
                nclp = currModel@nclp)
    if(correlation && rdf > 0) {
        ans$correlation <- (XtXinv * resvar)/outer(se, se)
        ans$symbolic.cor <- symbolic.cor
    }
    if(identical(object$call$algorithm, "port"))
	ans$message <- object$message
    class(ans) <- "summary.timp.nls"
    ans
}
print.summary.timp.nls <-
  function (x, digits = max(3, getOption("digits") - 3),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nFormula: ")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
    df <- x$df
    rdf <- df[2]
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 ...)
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    cat("Number of nonlinear parameters:", x$nnonlin, "\n")
    cat("Number of conditionally linear parameters:", x$nclp, "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Parameter Estimates:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {
		print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }

    .p.nls.convInfo(x, digits = digits)

    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep="")
    cat("\n")
    invisible(x)
}
print.summary.timp.nls.lm <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
            
{
  df <- x$df
  rdf <- df[2]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, ...)
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
  cat("Number of nonlinear parameters:", x$nnonlin, "\n")
  cat("Number of conditionally linear parameters:", x$nclp, "\n")
  cat("Number of iterations to termination:", x$niter, "\n")
  cat("Reason for termination:", x$stopmess, "\n")
  invisible(x)
}
.p.nls.convInfo <- function(x, digits)
{
    if(identical(x$call$algorithm, "port"))
	cat("\nAlgorithm \"port\", convergence message:",
	    x$message, "\n")
    else
	with(x$convInfo, {
	    cat("\nNumber of iterations",
		if(isConv) "to convergence:" else "till stop:", finIter,
		"\nAchieved convergence tolerance:",
                format(finTol, digits=digits),"\n")
	    if(!isConv)
		cat("Reason stopped:", stopMessage, "\n")
	})
    invisible()
}
print.summary.timp.optim <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
            
{
  df <- x$df
  rdf <- df[2]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, ...)
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
  cat("Number of nonlinear parameters:", x$nnonlin, "\n")
  cat("Number of conditionally linear parameters:", x$nclp, "\n")
  invisible(x)
}
print.timp.nls <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  
  cat("Nonlinear regression model\n")
  cat("  model: ", deparse(formula(x)), "\n")
  cat("   data: ", deparse(x$data), "\n")
  cat("Number of conditionally linear parameters:", x$nclp, "\n")
  cat("Number of nonlinear parameters:", length(x$m$getAllPars()), "\n")
  cat("nonlinear parameter estimates:\n")
  print(x$m$getAllPars(), digits = digits, ...)
  cat(" ", if(!is.null(x$weights) && diff(range(x$weights))) "weighted ",
      "residual sum-of-squares: ", format(x$m$deviance(), digits = digits),
      "\n", sep = '')
  .p.nls.convInfo(x, digits = digits)
  invisible(x)
}
print.timp.nls.lm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Nonlinear regression via the Levenberg-Marquardt algorithm\n")
    cat("Number of conditionally linear parameters:", x$nclp, "\n")
    cat("Number of nonlinear parameters:", length(x$par), "\n")
    cat("nonlinear parameter estimates:\n", toString(x$par), "\n")
    cat("residual sum-of-squares: ", format(x$deviance, digits = digits),
	"\n", sep = '')
    cat("reason terminated: ", x$message, "\n", sep='')
    invisible(x)
}
print.timp.optim <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Nonlinear regression via Nelder-Mead\n")
    cat("Number of conditionally linear parameters:", x$nclp, "\n")
    cat("Number of nonlinear parameters:", length(x$par), "\n")
    cat("nonlinear parameter estimates:\n", toString(x$par), "\n")
    invisible(x)
}
