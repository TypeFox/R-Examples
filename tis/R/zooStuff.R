jul.yearmon <- function(x, offset = 0, ...) jul(zoo::as.Date.yearmon(x, frac = offset))
ti.yearmon  <- function(x, ...) ti(jul(x), tif = "monthly")
jul.yearqtr <- function(x, offset = 0, ...) jul(zoo::as.Date.yearqtr(x, frac = offset))
ti.yearqtr  <- function(x, ...) ti(jul(x), tif = "quarterly")

as.tis.zoo <- function(x, ...){
  xindex <- attr(x, "index")  ## no need to import zoo:::index
  firstIndex <- xindex[1]
  if(is.ti(firstIndex)){
    attr(x, "index") <- NULL
    attr(x, "frequency") <- NULL
    return(tis(stripClass(stripClass(x, "zooreg"), "zoo"), start = firstIndex))
  }
  xti <- inferTi(xindex)
  zstart <- min(xti)
  zend   <- max(xti)
  if(is.matrix(x)){
    z <- tis(matrix(NA, nrow = zend - zstart + 1, ncol = ncol(x)), start = zstart)
    z[xti,] <- x[]
  }
  else {
    z <- tis(NA, start = zstart, end = zend)
    z[xti] <- x[]
  }
  z
}

