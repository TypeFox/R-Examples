print.weibull.frailty <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
   cat("\nFrailty distribution: Gamma")
   cat("\nFrailty variance:", round(x$coefficients$var.frailty, digits = digits))
   if (length(x$coefficients$betas)) {
        cat("\n\nCoefficients:\n")
        print(round(x$coefficients$betas, digits = digits))
   }
   cat("\nShape:", round(x$coefficients$shape, digits = digits))
   cat("\nScale:", round(x$coefficients$scale, digits = digits))
   cat("\n\nlog-Lik:", x$logLik, "\n\n")
   invisible(x)
}
