Vm.tfl <- function (obj, m, ...)
{
  if (! inherits(obj, "tfl")) stop("first argument must be object of class 'tfl'")
  m <- as.integer(m)
  if (length(m) != 1 || any(m < 0) ) stop("second argument must be single non-negative integer")

  if (attr(obj, "incomplete") && (m < attr(obj, "f.min") || m > attr(obj, "f.max"))) {
    NA
  }
  else {
    sum(obj$f == m)
  }
}
