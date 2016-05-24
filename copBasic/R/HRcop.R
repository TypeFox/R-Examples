"HRcop" <- function(u,v, para=NULL, ...) {
  if(para < 0) {
     warning("parameter must be 0 <= Theta < Infinity")
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

  if(para < .Machine$double.eps) return(P(u,v))
  d2 <- para/2; di <- 1/para
  x <- -log(u); y <- -log(v)
  opts <- options(warn=-1)
  cop <- exp(-x*pnorm(di + d2*log(x/y)) + -y*pnorm(di + d2*log(y/x)))
  cop[! is.finite(x)] <- 0; cop[! is.finite(y)] <- 0
  cop[x == 0] <- v # because v ---> 1
  cop[y == 0] <- u # because u ---> 1
  cop[is.nan(cop)] <- 1
  options(opts)
  return(cop)
}
