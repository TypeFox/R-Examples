## functions for prelimary R package wideLasso

expand.to.order.2 <- function(X, colsBinary, numBinary,
                              numNotBinary) {

  if(!is.matrix(X)) stop("X must be a matrix")
  if(storage.mode(X) != "double")
    stop("The storage mode of X must be double")

  return(.Call(expandToOrder2, X, as.integer(colsBinary),
               as.integer(numBinary), as.integer(numNotBinary)))

}

expand.to.order.3 <- function(X, X2, colsBinary, numBinary,
                              numNotBinary) {

  if(!is.matrix(X) || !is.matrix(X2))
    stop("X and X2 must be matrices")
  if(storage.mode(X) != "double" || storage.mode(X2) != "double")
    stop("The storage mode of X and X2 must be double")

  return(.Call(expandToOrder3, X, X2, as.integer(colsBinary),
               as.integer(numBinary), as.integer(numNotBinary)))

}

