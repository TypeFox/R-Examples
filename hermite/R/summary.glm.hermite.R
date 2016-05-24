summary.glm.hermite <- function (object, ...) 
{
  call <- attr(object, "Call")
  coef.p <- object$coefs
  covmat <- object$vcov
  var.cf <- diag(covmat)
  s.err <- sqrt(var.cf)
  s.err[length(coef.p)] <- NA
  tvalue <- coef.p/s.err
  tvalue[length(coef.p)] <- NA
  tvalue[length(coef.p)-1] <- object$w
  dn <- c("Estimate", "Std. Error")
  pvalue <- 2 * pnorm(-abs(tvalue))
  pvalue[length(coef.p)-1] <- object$pval
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "p-value"))
  aic <- 2*(length(coef.p)-1) - 2*object$loglik
  resid <- attr(object, "x") - object$fitted.values
  ans <- list(resid=resid, coefficients = coef.table, aic=aic, call=call)
  
  class(ans) <- "summary.glm.hermite"
  return(ans)
}