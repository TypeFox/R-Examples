Vm.spc <- function (obj, m, ...)
{
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  m <- as.integer(m)
  if ( any(is.na(m)) || any(m < 1) ) stop("second argument must be integer(s) >= 1")

  ## fill in 0's for empty frequency classes 
  idx <- match(m, obj$m)      
  Vm <- ifelse(is.na(idx), 0, obj$Vm[idx])

  ## unlisted classes in incomplete spectrum are set to NA
  m.max <- attr(obj, "m.max") 
  if (m.max > 0) {
    idx <- m > m.max
    if (any(idx)) Vm[idx] <- NA
  }

  Vm
}
