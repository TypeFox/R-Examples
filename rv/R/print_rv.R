# ========================================================================
# print.rv  -  print summary of a rv on the console
# ========================================================================

print.rv <- function(x, digits=rvpar("print.digits"), ...) {
  if (length(x)==0) {
    return(cat("rv(0)\n"))
  }
  print(summary(x, ...), digits=digits, ...)
}



