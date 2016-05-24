"asCOP" <- function(u,v, f=NULL, ...) {
  if(is.null(f)) {
    warning("'f' is a mandatory argument to be a function")
    return(NA)
  } else if(! is.function(f)) {
    warning("'f' must be a function")
    return(NA)
  }
  if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
    warning("length u = ",length(u), " and length v = ",length(v))
    warning("longer object length is not a multiple of shorter object length, no recycling")
    return(NA)
  }
  # The extra hassle of vectorization made here is to handle situations
  # in which nested integrals are used where uneven vectors can be passed
  if(length(u) == 1) {
     u <- rep(u, length(v))
  } else if(length(v) == 1) {
     v <- rep(v, length(u))
  }
  return(sapply(1:length(u), function(i) { f(u[i], v[i], ...) }))
}
