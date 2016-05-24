VVm.spc <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "spc")) stop("argument must be object of class 'spc'")
  m <- as.integer(m)
  if ( any(is.na(m)) || any(m < 1) ) stop("second argument must be integer(s) >= 1")
  if (!attr(obj, "hasVariances")) stop("no variance data available for frequency spectrum")
  if (!missing(N)) stop("argument 'N' not allowed for object of class 'spc'")

  ## fill in 0's for empty frequency classes 
  idx <- match(m, obj$m)      
  VVm <- ifelse(is.na(idx), 0, obj$VVm[idx])

  ## unlisted classes in incomplete spectrum are set to NA
  m.max <- attr(obj, "m.max") 
  if (m.max > 0) {
    idx <- m > m.max
    if (any(idx)) VVm[idx] <- NA
  }

  VVm
}
