print.summary.diagram <-
function(x, ...) {
  cat("Call: \n")
  print(x[["call"]])
  cat("\nNumber of features: \n")
  print(x[["n"]])
  cat("\nMax dimension: \n")
  print(x[["maxdimension"]])
  cat("\nScale: \n")
  print(x[["scale"]])
}