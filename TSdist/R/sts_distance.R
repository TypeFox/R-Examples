
# This function calculates the short time series (sts) distance measure.
STSDistance <- function(x, y, tx=NULL, ty=NULL) {
  
  # If no index is specified then evenly samples series are assumed.
  if (is.null(tx) & is.null(ty)) {
   tx <- c(1:length(x))
   ty <- tx
  }
  if (is.null(tx)) {
    tx <- ty
  }
  if (is.null(ty)) {
    ty <- tx
  }
  
  if (class(try(STSInitialCheck(x, y, tx, ty))) == "try-error") {
    return(NA)
  } else {

  # The STS distance is calculated.
  d <- sqrt(sum((diff(x) / diff(tx) - diff(y) / diff(ty)) ^ 2))
  return(d)
  }

}

#  This function checks for possible initial errors: 

STSInitialCheck <- function(x, y, tx, ty) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) <= 1 | length(y) <= 1) {    
    stop('The series must have a more than one point', call.=FALSE)
  }
  if (length(x) != length(y)) {    
    stop('Both series must have the same length', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  } 
  if (! missing(tx) & ! missing(ty)) {
    if (any(tx<=0) | any(ty<=0)) {      
      stop('The temporal indice must always be positive', call.=FALSE)
    }
    if (any(diff(tx) != diff(ty))) {      
      stop('The sampling rate must be equal in both series', call.=FALSE)
    }
    if (any(diff(tx) <= 0) | any(diff(ty) <= 0)) {      
      stop('The temporal index must be ascending.', call.=FALSE)
    }
    
    if ((length(tx) != length(x)) | (length(ty) != length(y))) {      
    stop('The length of the time indice must be equal to the length of the series', call.=FALSE)
  }
  } 
}