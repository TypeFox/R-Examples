summary.TML<-function(object, ...)
{
  coefs <- object$th1
  aliased <- is.na(coefs)
  if(any(aliased))
    coefs <- coefs[!aliased]
  p <- length(object$th1)
  n <- length(object$residuals)
  pi <- length(coefs)
  df <- c(p, n - p, pi)
  covmat <- object$COV
  if (is.null(covmat)) {
   cat("An estimate of the covariance matrix is required in object.\n",
       "Set cov='parametric' or cov='nonparametric'.\n") 
   return()}
  s.err <- sqrt(diag(covmat))[-(p + 1)]
  tvalue <- coefs/s.err
  pvalue <- 2 * pnorm(-abs(tvalue))
  coef.table <- cbind(coefs, s.err, tvalue, pvalue)
  colnames(coef.table) <- list("Estimate","Std.Error", "t-value", "Pr(>|t|)")
  ans <- c(object[c("call", "terms", "residuals", "fitted.values", "tn")],
    list(coefficients = coef.table, aliased = aliased, df = df, sigma = object$v1, 
    cutoff.values = c(object$tl, object$tu)), errors = object$errors)
  class(ans) <- "summary.TML"
  return(ans)
}


