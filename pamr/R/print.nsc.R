print.nsc <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  mat <- rbind(threshold = format(round(x$threshold, 3)), nonzero = 
               format(trunc(x$nonzero)), errors = x$errors)
  dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
  print(t(mat), quote = FALSE)
  invisible()
}
