N.spc <- function (obj, ...)
{
  if (! inherits(obj, "spc")) stop("argument must be object of class 'spc'")

  attr(obj, "N")
}
