nls.lm.control <- function(ftol = sqrt(.Machine$double.eps), ptol =
                           sqrt(.Machine$double.eps), gtol = 0, diag
                           = list(), epsfcn = 0, factor = 100, maxfev
                           = integer(), maxiter = 50, nprint = 0)
  list(ftol = ftol, ptol = ptol, gtol = gtol, diag = diag, epsfcn =
       epsfcn, factor = factor, maxfev = maxfev,
       maxiter = maxiter, nprint = nprint)
nls.lm <- function(par, lower=NULL, upper=NULL,fn, jac = NULL, 
       	  control = nls.lm.control(), ...)
{
  fn1  <- function(par) fn(par, ...)
  jac1 <- if (!is.null(jac))
    function(par) jac(par, ...)
  if(is.null(lower))
    lower <- rep(-Inf,length(par))
  if(is.null(upper))
    upper <- rep(Inf,length(par))
  if(length(lower) != length(par))
    stop("length(lower) must be equal to length(par)")
   if(length(upper) != length(par))
    stop("length(upper) must be equal to length(par)")
  ctrl <- nls.lm.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(length(ctrl[["maxfev"]]) == 0)
    ctrl[["maxfev"]] <- 100*(length(unlist(par)) + 1) 
  out <- .Call("nls_lm", par, lower, upper, fn1, jac1, ctrl, new.env(),
               PACKAGE = "minpack.lm")
  
  out$hessian <- matrix(out$hessian, nrow = length(unlist(par)))
  
  names(out$par)        <-
    rownames(out$hessian) <-
      colnames(out$hessian) <-
        names(out$diag)       <- names(par)
  out$deviance <- sum(out$fvec^2)
  class(out) <- "nls.lm"
  out
}
print.nls.lm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Nonlinear regression via the Levenberg-Marquardt algorithm\n")

    cat("parameter estimates:", toString(x$par), "\n")
    cat("residual sum-of-squares: ", format(x$deviance, digits = digits),
	"\n", sep = '')
    cat("reason terminated: ", x$message, "\n", sep='')
    invisible(x)
}
deviance.nls.lm <- function(object, ...) object$deviance
coef.nls.lm <- function(object, ...) unlist(object$par)
residuals.nls.lm <- function(object, ...) object$fvec
df.residual.nls.lm <- function(object, ...)
  length(resid(object)) - length(coef(object))
summary.nls.lm <- function (object, ...)
{
    param <- coef(object)
    pnames <- names(param)
    ibb <- chol(object$hessian)
    ih <- chol2inv(ibb)
    p <- length(param)
    rdf <- length(object$fvec) - p 
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
                coefficients = param)
    class(ans) <- "summary.nls.lm"
    ans
}
print.summary.nls.lm <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
            
{
  df <- x$df
  rdf <- df[2]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, ...)
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
  cat("Number of iterations to termination:", x$niter, "\n")
  cat("Reason for termination:", x$stopmess, "\n")
  invisible(x)
}

vcov.nls.lm <- function(object,...) {
  object$deviance/(length(object$fvec)-length(object$par))*solve(object$hessian)
}

confint.nls.lm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(object)
  if (missing(parm)) parm <- seq_along(cc)
  levs <- c((1-level)/2,0.5+level/2)
  dfval <- (length(object$fvec)-length(object$par))
  tdist <- qt(levs[2],dfval)
  m1 <- outer(sqrt(diag(vcov(object))),c(-1,1))*tdist
  m2 <- sweep(m1,cc,MARGIN=1,"+")
  colnames(m2) <- paste(100*levs,"%",sep="")
  m2[parm,]
}

