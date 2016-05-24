V.tfl <- function (obj, ...)
{
  if (! inherits(obj, "tfl")) stop("argument must be object of class 'tfl'")

  attr(obj, "V")
}
