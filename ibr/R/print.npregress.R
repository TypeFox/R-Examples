print.npregress <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nDf:",format(round(x$df,digits-1)),"; bandwidth:", format(round(x$bandwidth,digits)), "\n")
  cat("Chosen by:", x$call$criterion, "\n")
}
