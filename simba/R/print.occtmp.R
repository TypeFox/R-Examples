"print.occtmp" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\n")
  cat("Change in species occurrence", "\n")
  print(x$bac, digits = 3)
  cat("\n\n")
  cat("Further Statistics:", "\n")
  print(x$stats, digits = 3)
  cat("\n\n")
  invisible(x)
}