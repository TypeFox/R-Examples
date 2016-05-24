Vm.vgc <- function (obj, m, ...)
{
  if (! inherits(obj, "vgc")) stop("first argument must be object of class 'vgc'")
  m <- as.integer(m)
  if ( length(m) != 1 || any(is.na(m)) || any(m < 1) ) stop("second argument must be integer >= 1")

  varname <- paste("V", m, sep="")
  if (m > attr(obj, "m.max"))
    stop("spectrum element '",varname,"' not included in vocabulary growth curve")
  obj[[ varname ]]
}
