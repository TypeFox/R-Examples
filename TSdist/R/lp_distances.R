
# Wrapper function for all lp distances.
LPDistance <- function(x, y, method="euclidean", ...) {
    
  method <- match.arg(method, c("euclidean", "minkowski", "infnorm", "manhattan"))
  
  if (method == "euclidean") {
    d <- EuclideanDistance(x, y)
  }
  if (method == "minkowski") {
    d <- MinkowskiDistance(x, y, ...)
  }
  if (method == "infnorm") {
    d <- InfNormDistance(x, y)
  }
  if (method == "manhattan") {
    d <- ManhattanDistance(x, y)
  }
  
  return(d)
}

# This function checks for initial errors.
LPInitialCheck <- function(x, y, p) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric.', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) < 1 | length(y) < 1) {
    stop('The series must have at least one point.', call.=FALSE)
  }
  if (length(x) != length(y)) {
    stop('Both series must have the same length.', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series.', call.=FALSE)
  } 
  if (! missing(p)) {
    
    if (round(p) != p) {
      stop('p must be an integer value.', call.=FALSE)
    }
    if (p <= 0) {
      stop('p must be positive.', call.=FALSE)
    }
  } 
}