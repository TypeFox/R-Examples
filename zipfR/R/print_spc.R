print.spc <- function (x, all=FALSE, ...)
{
  if (! inherits(x, "spc")) stop("argument must be object of class 'spc'")

  complete <- attr(x, "m.max") == 0
  if (!all && nrow(x) > 10) {
    x <- x[1:10,]
    complete <- FALSE
  }

  print.data.frame(x)
  if (!complete) cat("\t...\n")
  cat("\n")
  
  info <- data.frame(N=attr(x,"N"), V=attr(x,"V"))
  if (attr(x, "hasVariances")) info$VV <- attr(x, "VV")
  rownames(info) <- ""
  print(info)

  invisible(NULL)
}
