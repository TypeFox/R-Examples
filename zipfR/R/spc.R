spc <- function (Vm, m=1:length(Vm), VVm=NULL, N=NA, V=NA, VV=NA,
                 m.max=0, expected=!missing(VVm))
{
  if (length(m) != length(Vm)) stop("'m' and 'Vm' must have the same length")
  if (missing(N) != missing(V)) stop("'N' and 'V' must always be specified together")
  if (missing(m.max) && !missing(N)) {
    if (sum(Vm) != V || sum(m*Vm) != N) m.max <- max(m) # spectrum seems to be incomplete
  }
  incomplete <- (m.max > 0)
  if (incomplete && missing(N)) stop("'N' and 'V' must be specified for incomplete spectrum")
  variances <- !missing(VVm)
  if (variances && !expected) stop("'VVm' may only be specified for expected frequency spectrum")
  if (variances && length(VVm) != length(Vm)) stop("'VVm' must have same length as 'm' and 'Vm'")
  if (variances && missing(VV)) warning("if variances are present, 'VV' should also be specified")
  if (! is.integer(m)) {
    if (any(m != floor(m))) stop("class sizes 'm' must be integer values")
    m <- floor(m)
  }

  m <- as.double(m)         # make sure there are no integer overflows
  Vm <- as.double(Vm)
  if (variances) VVm <- as.double(VVm)
  
  # remove empty frequency classes for compact storage
  idx <- (Vm != 0)
  if (any(!idx)) {
    m <- m[idx]
    Vm <- Vm[idx]
    if (variances) VVm <- VVm[idx]
  }

  # calculate N anv V automatically (for complete spectrum)
  if (missing(N)) {
    N <- sum(m * Vm)
    V <- sum(Vm)
  } else {
    N <- as.double(N)
    V <- as.double(V)
  }

  # remove any frequency classes above m.max from the spectrum
  if (incomplete && max(m) > m.max) {
    idx <- (m <= m.max)
    m <- m[idx]
    Vm <- Vm[idx]
    if (variances) VVm <- VVm[idx]
  }

  # consistency check for incomplete frequency spectrum
  if (incomplete) {
    miss.N <- N - sum(m * Vm)
    miss.V <- V - sum(Vm)
    if (miss.V < -1e-12 * V) stop("inconsistent data (V=",V," < sum(Vm)=",sum(Vm),")")
    # tolerant check for error condition miss.V * (m.max+1) > miss.N (because of rounding errors in expected frequency spectra)
    if (miss.V * (m.max+1) - miss.N > 1e-12 * N) {
      stop("inconsistent data (N=",N," should be at least ", sum(m * Vm) + miss.V * (m.max+1),")")
    }
  }

  # make sure that spectrum elements are in increasing order
  idx <- order(m)
  if (any(diff(m) != 1)) {
    m <- m[idx]
    Vm <- Vm[idx]
    if (variances) VVm <- VVm[idx]
  }
  
  spc <- data.frame(m=m, Vm=Vm)
  if (variances) spc$VVm <- VVm
  attr(spc, "m.max") <- m.max
  attr(spc, "N") <- N
  attr(spc, "V") <- V
  attr(spc, "expected") <- expected
  attr(spc, "hasVariances") <- variances
  if (variances) attr(spc, "VV") <- as.double(VV)
  
  class(spc) <- c("spc", class(spc))
  spc
}
