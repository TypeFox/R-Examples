print.summary.maxPersistence <-
function(x, ...) {
  cat("Call: \n")
  print(x[["call"]])
  cat("\nThe number of significant features is maximized by \n")
  print(x[["maxNumParam"]])
  cat("\nThe total significant persistence is maximized by \n")
  print(x[["maxSigParam"]])
}