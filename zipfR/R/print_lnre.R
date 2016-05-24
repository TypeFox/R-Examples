print.lnre <- function (x, ...)
{
  if (! inherits(x, "lnre")) stop("argument must belong to a subclass of 'lnre'")

  summary.lnre(x)
}
