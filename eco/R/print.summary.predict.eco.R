print.summary.predict.eco <- function(x, digits=max(3, getOption("digits")
                                           -3), ...) {
  cat("\nOut-of-sample Prediction:\n")
  print(x$W.table, digits=digits, na.print="NA",...)

  cat("\nNumber of Monte Carlo Draws:", x$n.draws)
  
  cat("\n")
  invisible(x)
}
