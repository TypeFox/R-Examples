VV.vgc <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "vgc")) stop("argument must be object of class 'vgc'")
  if (!attr(obj, "hasVariances")) stop("no variance data available for vocabulary growth curve")
  if (!missing(N)) stop("argument 'N' not allowed for object of class 'vgc'")
  
  obj$VV
}
