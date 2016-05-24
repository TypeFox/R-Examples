print.summary.oneway.rfit <- function (x, digits = max(5, .Options$digits - 2),
  ...) {
  cat("\nMultiple Comparisons\n")
  cat(paste("Method Used ",x$method),"\n\n")
  print(cbind(x$table[,1:2],round(x$table[,3:6],digits)))
}

