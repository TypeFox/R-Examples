vgc <- function (N, V, Vm=NULL, VV=NULL, VVm=NULL, expected=FALSE, check=TRUE)
{
  if (length(N) != length(V)) stop("'N' and 'V' must have the same length")
  if (any(N < 0)) stop("sample sizes 'N' must be non-negative integers")
  if (check) {
    if (any(diff(N) < 0)) stop("sample sizes 'N' must be increasing")
    if (any(V < 0)) stop("vocabulary sizes 'V' must be non-negative")
    if (any(diff(V) < 0)) stop("inconsistent decrease in vocabulary size 'V'")
  }

  N <- as.double(N)             # make sure to avoid integer overflows
  V <- as.double(V)

  if (!missing(Vm)) {
    if (!is.list(Vm)) Vm <- as.list(Vm) # allow V1 to be specified as plain vector
    for (v.m in Vm) {
      if (!(is.numeric(v.m))) stop("elements of 'Vm' must be numeric vectors")
      if (length(v.m) != length(N)) stop("vectors in 'Vm' must have the same length as 'N' and 'V'")
      if (check) {
        if (any(v.m < 0)) stop("spectrum elements in 'Vm' must be non-negative")
      }
    }
    Vm <- lapply(Vm, as.double) # make sure to avoid integer overflows
  }
  m.max <- length(Vm)               # m.max = 0 if Vm is not specified
  ## -- this limitation seems arbitrary, so it has been disabled
  ## if (m.max > 9) stop("at most 9 spectrum elements allowed in 'Vm'")

  variances <- !missing(VV)
  if (variances && !is.list(VVm)) VVm <- list(VVm) # same as above
  m.max.var <- length(VVm)              # spectrum elements for which variances are specified
  if (variances && (m.max.var != m.max))
    stop("variances 'VVm' must be specified for the same spectrum elements as in 'Vm'")
  if (!variances && m.max.var > 0)
    stop("variance 'VV' missing (but variances 'VVm' are specified)")
  if (variances) {
    for (vv.m in VVm) {
      if (!is.numeric(vv.m)) stop("elements of 'VVm' must be numeric vectors")
      if (length(vv.m) != length(N)) stop("vectors in 'VVm' must have the same length as 'N' and 'V'")
      if (check) {
        if (any(vv.m < 0)) stop("variances in 'VVm' must be non-negative")
      }
    }
    VVm <- lapply(VVm, as.double)
  }
  
  vgc <- data.frame(N=N, V=V)
  if (variances) vgc$VV <- VV
  if (m.max > 0) {
    for (i in 1:m.max) {
      vgc[[ paste("V", i, sep="") ]] <- Vm[[i]]
      if (variances) vgc[[ paste("VV", i, sep="") ]] <- VVm[[i]]
    }
  }
    
  attr(vgc, "m.max") <- m.max
  attr(vgc, "expected") <- expected || variances
  attr(vgc, "hasVariances") <- variances
  
  class(vgc) <- c("vgc", class(vgc))
  vgc
}
