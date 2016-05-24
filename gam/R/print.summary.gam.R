"print.summary.gam" <-
  function(x,  digits = max(3, getOption("digits") - 3), quote = TRUE, prefix = "", ...)
{
  cat("\nCall: ")
  dput(x$call)
  dresid <- x$deviance.resid
  n <- length(dresid)
  rdf <- x$df[2]
  if(rdf > 5) {
    cat("Deviance Residuals:\n")
    rq <- quantile(as.vector(dresid))
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  }
  else if(rdf > 0) {
    cat("Deviance Residuals:\n")
    print(dresid, digits = digits)
  }
  cat(paste("\n(Dispersion Parameter for ", names(x$dispersion), 
            " family taken to be ", format(round(x$dispersion, digits)),
            ")\n",sep=""))
  int <- attr(x$terms, "intercept")
  cat("\n    Null Deviance:", format(round(x$null.deviance, digits)),
      "on", n - int, "degrees of freedom")
  cat("\nResidual Deviance:", format(round(x$deviance, digits)), "on",
      format(round(rdf, digits)), "degrees of freedom")
  cat("\nAIC:", format(round(x$aic, digits)),"\n")
  if(!is.null(x$na.action))
    cat(naprint(x$na.action), "\n")
  cat("\nNumber of Local Scoring Iterations:", format(trunc(x$iter)),
      "\n")
  aod=x$parametric.anova
  cat("\n")
  if(!is.null(aod)) print(aod)
  aod=x$anova
  cat("\n")
  if(!is.null(aod)) print(aod)
}
