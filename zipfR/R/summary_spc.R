summary.spc <- function (object, ...)
{
  if (! inherits(object, "spc")) stop("argument must be object of class 'spc'")
  m.max <- attr(object, "m.max")
  K <- 8                       # show first K spectrum elements (V_1 ... V_K)
  
  cat("zipfR object for ")
  if (attr(object, "expected")) cat("expected ")
  cat("frequency spectrum")
  if (m.max > 0) cat(", incomplete (m <= ", m.max, ")", sep="")
  if (attr(object, "hasVariances")) cat(", with variances")
  cat("\n")

  cat("Sample size:     N  =", attr(object, "N"), "\n")
  cat("Vocabulary size: V  =", attr(object, "V"), "\n")
  cat("Class sizes:     Vm = ")
  if (m.max > 0 && m.max < K) {
    Vm <- Vm(object, m=1:m.max)
  }
  else {
    Vm <- Vm(object, m=1:K)
  }
  cat(Vm)
  if (m.max > 0 || max(object$m) > K) cat(" ...")
  cat("\n")

  invisible(NULL)
}
