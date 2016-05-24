densitynormalmulti <- function(x, bw=NULL, V=NULL) {
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(bw)) {
    bw <- (4/(2*p+1))^(1/(p+4))*n^(-1/(p+4))
  }
  if (is.null(V)) {
    V <- cov(x)
  }
  detV <- det(V)
  V <- solve(V)/(bw^2)
  z <- .Fortran("dennormu",
    as.double(x),
    as.integer(n),
    as.integer(p),
    as.double(V),
    den=double(n),
    PACKAGE="wle"            
  )
  den <- z$den/(n*sqrt(detV)*(2*pi)^(p/2)*bw^p)
  return(den)
}
