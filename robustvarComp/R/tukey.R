dtukeychi <- function(x, cc) {
  n <- length(x)
  z <- 0
  z <- .Fortran("dtukeych",
    as.double(x),
    as.integer(n),
    as.double(cc),
    z=as.double(z),
    PACKAGE="robustvarComp")$z
  return(z)  
}
