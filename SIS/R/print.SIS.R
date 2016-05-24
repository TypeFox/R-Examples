print.SIS <- function(x, digits=max(3, getOption("digits") - 3), ...){
   cat("\nCall: ", deparse(x$call), "\n\n")
   cat("\n$ix", "\n\n")
   print(signif(x$ix, digits))
   cat("\n$coef.est", "\n\n")
   print(signif(x$coef.est, digits))
}
