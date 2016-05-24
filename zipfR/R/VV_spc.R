VV.spc <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "spc")) stop("argument must be object of class 'spc'")
  if (!attr(obj, "hasVariances")) stop("no variance data available for frequency spectrum")
  if (!missing(N)) stop("argument 'N' not allowed for object of class 'spc'")
  
  attr(obj, "VV")
}
