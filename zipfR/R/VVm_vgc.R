VVm.vgc <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "vgc")) stop("argument must be object of class 'vgc'")
  m <- as.integer(m)
  if ( length(m) != 1 || any(is.na(m)) || any(m < 1) ) stop("second argument must be integer >= 1")
  if (!attr(obj, "hasVariances")) stop("no variance data available for vocabulary growth curve")
  if (!missing(N)) stop("argument 'N' not allowed for object of class 'vgc'")

  varname <- paste("VV", m, sep="")
  if (m > attr(obj, "m.max"))
    stop("spectrum variances '",varname,"' not included in vocabulary growth curve")
  obj[[ varname ]]
}
