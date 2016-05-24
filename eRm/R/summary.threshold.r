summary.threshold <- function(object,...)
{
#object of class "threshold"

  coef.table <- cbind(round(object$threshpar,5),round(object$se.thresh,5),round(confint(object),5))
  dimnames(coef.table) <- list(names(object$threshpar),c("Estimate","Std. Err.",colnames(confint(object))))
  cat("\n")
  print(coef.table)
  cat("\n")
}