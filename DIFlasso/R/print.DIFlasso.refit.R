print.DIFlasso.refit <-
function(x, ...){
  cat("Number of (valid) persons: P =",x$P,"\n")
  cat("Number of items: I =",x$I,"\n")
  cat("DIF Items:",x$dif.items,"\n")
  cat("\n")
  cat("Matrix of estimated item-specific coefficients (after refit):\n")
  print.default(x$gamma, ...)
  invisible(x)
}
