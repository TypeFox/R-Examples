print.flexCPH <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   print(lapply(x$coefficients, round, digits = digits))
   cat("\nlog-Lik:", x$logLik, "\n\n")
   invisible(x)
}
