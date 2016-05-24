tfl <- function  (f, k=1:length(f), type=NULL, f.min=min(f), f.max=max(f),
                  incomplete=!(missing(f.min) && missing(f.max)), N=NA, V=NA,
                  delete.zeros=FALSE)
{
  if (! (is.numeric(f) && all(f >= 0))) stop("'f' must be non-negative integer vector")
  if (length(f) == 0 && missing(k)) k <- integer(0)   # since 1:0 = c(1,0)
  if (length(k) != length(f)) stop("'k' and 'f' must have the same length")
  hasTypes <- !is.null(type)
  if (hasTypes && length(type) != length(k)) stop("'type' must have the same length as 'k' and 'f'")

  f <- as.double(f)         # make sure there are no integer overflows
  
  if (delete.zeros) {
    idx <- f == 0
    if (any(idx)) {
      f <- f[!idx]
      k <- k[!idx]
      if (hasTypes) type <- type[!idx]
    }
  }

  if (missing(N) != missing(V)) stop("'N' and 'V' must always be specified together")
  if (missing(N)) {
    if (incomplete) stop("'N' and 'V' must be specified for incomplete spectrum")
    N <- sum(f)
    V <- sum(f > 0)
  } else {
    N <- as.double(N)
    V <- as.double(V)
  } 
  
  ## limit frequencies to specified range
  if (missing(f.min) & length(f) == 0) f.min <- 0 # avoid warnings for empty TFLs 
  if (missing(f.max) & length(f) == 0) f.max <- 0
  idx <- f.min <= f & f <= f.max
  if (any(!idx)) {
    k <- k[idx]
    f <- f[idx]
    if (hasTypes) type <- type[idx]
    incomplete <- TRUE # now tfl is known to be incomplete, even if input data wasn't
  }

  ## ensure that rows are sorted properly (IDs 'k' must be in ascending order)
  idx <- order(k)
  if (any(diff(idx) != 1)) {
    k <- k[idx]
    f <- f[idx]
    if (hasTypes) type <- type[idx]
  }
  
  tfl <- data.frame(k=k, f=f)
  if (hasTypes) tfl$type <- type
  attr(tfl, "N") <- N
  attr(tfl, "V") <- V
  attr(tfl, "f.min") <- f.min
  attr(tfl, "f.max") <- f.max
  attr(tfl, "incomplete") <- incomplete
  attr(tfl, "hasTypes") <- hasTypes

  class(tfl) <- c("tfl", class(tfl))
  tfl
}
