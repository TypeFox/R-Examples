optimChi <- function(x, cc) {
  n <- length(x)
  y <- 0
  z <- .Fortran("optimchi",
    as.double(x),
    as.double(cc),
    y=as.double(y),
    PACKAGE="robustvarComp")$y
  return(z)
}

doptimchi <- function(x, cc) {
  n <- length(x)
  z <- 0
  z <- .Fortran("doptimch",
    as.double(x),
    as.integer(n),
    as.double(cc),
    z=as.double(z),
    PACKAGE="robustvarComp")$z
  return(z)
}
