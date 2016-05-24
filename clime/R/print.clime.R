print.clime <- function(x,digits = max(3, getOption("digits") - 3), ... ) {
  cat("\n clime options summary: \n")
  cat(" perturb=", signif(x$perturb, digits), "\t standardize=", x$standardize)
  cat("\n lambdas used:\n")
  print(signif(x$lambda,digits))
}
