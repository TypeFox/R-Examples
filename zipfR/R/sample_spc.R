sample.spc <- function (obj, N, force.list=FALSE) {
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  if (attr(obj, "m.max") > 0) stop("incomplete frequency spectra are not supported")
  if (attr(obj, "expected")) stop("expected frequency spectra are not supported")
  if (!is.integer(obj$Vm) && any(obj$Vm != floor(obj$Vm)))
    stop("all spectrum elements V_m must be integer values")
  if (! (is.numeric(N) && all(N >= 0))) stop("N must be vector of non-negative integers")
  if (any(N > N(obj))) stop("can't sample ", max(N), " tokens from 'spc' object with N=", N(obj))

  result <- sample.tfl(spc2tfl(obj), N, force.list=TRUE)
  result <- lapply(result, tfl2spc)
  
  if (length(N) == 1 && !force.list) {
    result[[ 1 ]]
  }
  else {
    result
  }
}
