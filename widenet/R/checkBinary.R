## wrapper for C function which checks columns of X to see
## which are binary

check.binary <- function(X) {

  if(!is.matrix(X)) stop("X must be a matrix")
  if(storage.mode(X) != "double")
    stop("The storage mode of X must be double")

  return(.Call(checkColsBinary, X))
}
