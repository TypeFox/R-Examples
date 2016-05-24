summary.vgc <- function (object, ...)
{
  if (! inherits(object, "vgc")) stop("argument must be object of class 'vgc'")

  cat("zipfR object for ")
  if (attr(object, "expected")) cat("expected ")
  cat("vocabulary growth curve")
  if (attr(object, "hasVariances")) cat(", with variances")
  cat("\n")

  m.max <- attr(object, "m.max")
  cat(nrow(object), "samples for N =", min(object$N), "...", max(object$N), "\n")
  if (m.max > 0) cat("Spectrum elements included up to m =", m.max, "\n")

  invisible(NULL)
}
