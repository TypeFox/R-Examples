print.vgc <- function (x, all=FALSE, ...)
{
  if (! inherits(x, "vgc")) stop("argument must be object of class 'vgc'")

  n <- nrow(x)
  if (!all && n > 25) {
    idx <- sort(sample(n, 25))
    print.data.frame(x[idx,])
    cat("\t(random subset of 25 entries shown)\n")
  }
  else {
    print.data.frame(x)
  }

  invisible(NULL)
}
