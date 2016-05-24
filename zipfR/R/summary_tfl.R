summary.tfl <- function (object, ...)
{
  if (! inherits(object, "tfl")) stop("argument must be object of class 'tfl'")

  incomplete <- attr(object, "incomplete")
  f.min <- attr(object, "f.min")
  f.max <- attr(object, "f.max")
  f <- object$f
  
  cat("zipfR object for ")
  if (incomplete) cat("incomplete ")
  cat("frequency spectrum")
  if (incomplete) cat(", restricted to range", f.min, "...", f.max)
  cat("\n")

  cat("Sample size:     N  =", attr(object, "N"), "\n")
  cat("Vocabulary size: V  =", attr(object, "V"), "\n")
  if (!incomplete) {
    cat("Range of freq's: f  =", f.min, "...", f.max, "\n")
    cat("Mean / median:   mu =", mean(f), ",  M =", median(f), "\n")
    cat("Hapaxes etc.:    V1 =", sum(f==1), ",  V2 =", sum(f==2), "\n")
  }
  if (attr(object, "hasTypes")) {
    idx <- order(object$type)
    all.shown <- !(length(idx) > 6)
    if (!all.shown) idx <- idx[1:6]
    cat("Types:  ", as.character(object$type[idx]))
    if (!all.shown) cat(" ...")
    cat("\n")
  }

  invisible(NULL)
}
