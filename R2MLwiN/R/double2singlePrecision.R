#' Converts numerical values from double precision to single precision.
#' 
#' This function converts numeric column(s) of a data frame object, matrix or
#' vector from double precision to single precision, e.g. to avoid a warning
#' from MLwiN which currently only stores data in single precision.
#' 
#' @param x A \code{\link[base]{data.frame}} object, matrix or vector to be converted. Column(s)
#' of these objects will be ignored during conversion if they are not
#' numeric.
#' 
#' @return An object of numerical values in single precision will be returned.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#' 
#' @export
double2singlePrecision <- function(x) {
  ## A wrapped function is used to change double to single
  if (is.matrix(x) || is.data.frame(x)) {
    flag <- 0
    if (is.data.frame(x)) 
      flag <- 1
    x <- apply(x, 2, function(y) if (is.double(y)) 
      .changePrecision(y, size = 4) else y)
    if (flag == 1) 
      x <- as.data.frame(x)
  } else {
    if (is.vector(x) && is.double(x)) {
      x <- .changePrecision(x, size = 4)
    }
  }
  x
}

### function .changePrecision
##  converts double values to double values in a given precision
##  (only correctly working for cut a higher precision to a lower one; e.g.
##  IEEE 754 double precision to IEEE 754 single precision)
.changePrecision <- function(x, size) {
  ## Sebastian Gibb
  
  # create a raw object to avoid direct file access
  virtualCon <- raw()
  # write binary data to raw object and change (mostly cut) precision to size size==4 # 32bit, single precision
  # size==8 # 64bit, double precision
  virtualCon <- writeBin(object = x, con = virtualCon, size = size)
  # re-read data
  x <- readBin(con = virtualCon, what = double(), size = size, n = length(x))
  return(x)
} 
