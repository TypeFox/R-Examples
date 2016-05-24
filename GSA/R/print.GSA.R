print.GSA <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  mat1=x$fdr.lo
print("")
print("Negative")
  print(mat1, quote = FALSE)
  mat2=x$fdr.hi
print("")
print("")
print("Positive")
  print(mat2, quote = FALSE)

  invisible()
}

