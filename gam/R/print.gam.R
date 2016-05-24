"print.gam" <-
  function(x, digits = 5, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  n <- x$df.null
  if(is.null(df.resid <- x$df.resid))
    df.resid <- n - sum(!is.na(x$coef)) - sum(x$nl.df)
  cat("\nDegrees of Freedom:", n, "total;", format(round(df.resid, digits
                                                         )), "Residual\n")
  if(!is.null(x$na.action))
    cat(naprint(x$na.action), "\n")
  cat("Residual Deviance:", format(round(x$deviance, digits)), "\n")
  invisible(x)
}
