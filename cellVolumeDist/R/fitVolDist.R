`fitVolDist` <-
function(vol, freq, r=100,sigma_r=44, t=40,sigma_t=.3*t,
                       maxiter=100, nprint=1, alg="leastsq"){

  if(alg=="leastsq") {
    nls.lm.out <- nls.lm(par=list(r=r,sigma_r=sigma_r),
                         fn=resFun, volFreq = freq,
                         t=t, sigma_t=.3*t, V=vol,
                         control=nls.lm.control(maxiter=maxiter, nprint=1)) 
    fv <- fitFun(nls.lm.out$par, sigma_t=sigma_t, t=t, V=vol)  
    out <- list(t=t, sigma_t = sigma_t, fitted=fv,
                fit.res = nls.lm.out, data = freq, 
                summary.fit.res = summary(nls.lm.out), alg="leastsq")    
  }
  else if(alg=="chisq") {
    optim.out <- optim(par=c(r, sigma_r),
                       fn=resFunChiSq, volFreq = freq,
                       t=t, sigma_t=.3*t, V=vol,
                       skel = list(r=r,sigma_r=sigma_r),
                       method = "BFGS",
                       control=list(maxit=maxiter, trace=1))
    fv <- fitFun(relist(optim.out$par, skeleton = list(r=r,sigma_r=sigma_r)),
                 sigma_t=sigma_t, t=t, V=vol)
    out <- list(t=t, sigma_t = sigma_t, fitted=fv, data = freq, 
                fit.res = optim.out, alg="chisq")
  }
  else{
    warning("Argument `alg' is invalid; specify either `leastsq' or `chisq'")
    return()
  }
  class(out) <- "fitVolDist"
  out
}
`print.fitVolDist` <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Volume distrubution model fit under least squares criteria\n")

  cat("Fixed model parameters:\n")
  cat("t=",x$t,"sigma_t=",x$sigma_t,"\n")
  ## following lines are from summary.nls.lm
  if(x$alg == "leastsq") {
    object <- x$fit.res 
    param <- coef(object)
    pnames <- names(param)
    p <- length(param)
    rdf <- length(object$fvec) - p 
    resvar <- deviance(object) / rdf
    se <- sqrt(diag(x$summary.fit.res$cov.unscaled) * resvar)
    names(se) <- pnames
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    cat("\nFree model parameters:\n")
    printCoefmat(param, digits = digits, ...)
    cat("Number of iterations to termination:", x$fit.res$niter, "\n")
    cat("Residual sum of squares:", x$fit.res$deviance, "\n")
  }
  if(x$alg == "chisq") {
    cat("Parameter estimates:\n")
    pp <- coef(x) 
    cat("r:", pp[1], "\n")
    cat("sigma_r",  pp[2], "\n")
  }
  invisible(x)
}
fitted.fitVolDist <- function(object, ...)  object$fitted
residuals.fitVolDist <- function(object, ...)  {
  if(object$alg == "leastsq") 
    x <- object$fit.res$fvec
  else if(object$alg == "chisq")
    x <- object$data - object$fv
  x
}
coef.fitVolDist <- function(object, ...) unlist(object$fit.res$par)
deviance.fitVolDist <- function(object, ...) {
   if(object$alg == "leastsq") 
    x <- object$fit.res$deviance  
  else if(object$alg == "chisq")
    x <- sum((object$data - object$fv)^2)
   x
}
  
vcov.fitVolDist <- function(object, ...)
{
  if(object$alg == "leastsq") {
     sm <- object$summary.fit.res
     x <- sm$cov.unscaled * sm$sigma^2
  }
  else if(object$alg == "chisq") {
    x <- NULL
    warning("Covariance matrix is not available")
  }
  x
  
 
}
