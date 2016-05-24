#' @export
#' @importFrom stats printCoefmat pnorm
#' @importFrom methods setClass setMethod

summary.aeefit <- function(object, digits = 2, dig.tst = 4, ...) {
  #if(is.null(digits)) digits <- options()$digits
  #else options(digits = digits)
  cat("\nCall:\n");   print(object$formula)
  cat("\nCoefficients:\n");  
  z <- object$beta/object$betaSE
  pval <- (pnorm(-abs(z)))*2
  cmat <- cbind(object$beta,exp(object$beta),object$betaSE,as.vector(z),as.vector(pval))              
  colnames(cmat) <- c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")
  printCoefmat(cmat, digits, dig.tst,
               signif.stars = TRUE,
               signif.legend = TRUE)
}
