print.fittedloop <- function (x,...) {
  cat("Call:\n")
  print(x$call)
  cat("Estimates and Delta Method Standard Errors:\n")
  print(cbind("Estimates"=x$values,"Std.Errors"=x$Std.Errors))
}
