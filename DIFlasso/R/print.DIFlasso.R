print.DIFlasso <-
function(x, ...){
  I <- nrow(x$beta)
  dif.mat <- x$dif.mat
  dif.items <- x$dif.items
  cat("Number of (valid) persons: P =",x$P,"\n")
  cat("Number of items: I =",I,"\n")
  cat("DIF Items:",dif.items,"\n")
  cat("\n")
  cat("Matrix of estimated item-specific coefficients:\n")
  print.default(dif.mat, ...)
  invisible(x)
}
