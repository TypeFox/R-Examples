print.tfl <- function (x, all=FALSE, ...)
{
  if (! inherits(x, "tfl")) stop("argument must be object of class 'tfl'")

  idx <- order(x$f, decreasing=TRUE)
  all.shown <- TRUE
  if (!all && length(idx) > 20) {
    idx <- idx[1:20]
    all.shown <- FALSE
  }

  print.data.frame(x[idx,])
  if (!all.shown) cat("\t...\n")
  cat("\n")
  
  info <- data.frame(N=attr(x,"N"), V=attr(x,"V"))
  if (attr(x,"incomplete")) {
    info$f.min <- attr(x, "f.min")
    info$f.max <- attr(x, "f.max")
  }
  rownames(info) <- ""
  print(info)

  invisible(NULL)
}
