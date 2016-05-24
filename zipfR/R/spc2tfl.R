spc2tfl <- function (spc)
{
  if (! inherits(spc, "spc")) stop("argument must be object of class 'spc'")
  if (attr(spc, "m.max") > 0) stop("incomplete frequency spectra are not supported")
  if (attr(spc, "expected")) stop("expected frequency spectra are not supported")
  if (!is.integer(spc$m)) stop("frequency classes in 'spc' must be integer values")
  if (!is.integer(spc$Vm) && any(spc$Vm != floor(spc$Vm)))
    stop("class sizes in 'spc' must be integer values")
  
  x <- list(values=rev(spc$m), lengths=rev(as.integer(spc$Vm)))
  class(x) <- "rle"
  tfl(inverse.rle(x))
}
