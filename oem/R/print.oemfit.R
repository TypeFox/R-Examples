print.oemfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call))
  cat("\nPenalty:", x$penalty, "\n\n")
  print(cbind(Df = x$df, sumSquare = signif(x$sumSquare, digits),
              Lambda = signif(x$lambda, digits)))
}
