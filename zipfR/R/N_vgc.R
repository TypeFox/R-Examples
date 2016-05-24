N.vgc <- function (obj, ...)
{
  if (! inherits(obj, "vgc")) stop("argument must be object of class 'vgc'")

  obj$N
}
