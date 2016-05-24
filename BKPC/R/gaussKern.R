gaussKern <-
function(x, newdata = x, theta = NULL){
  nr <- dim(newdata)[1]
  nc <- dim(x)[1]
  nf <- dim(x)[2]
  if (nf!= dim(newdata)[2])stop("error: newdata does not have the same number of features as data")
  dot <- .C("rkern", as.double(t(newdata)), as.double(t(x)), K = as.double(matrix(0, nr * nc, 1) ),  as.integer(nr), as.integer(nc), as.integer(nf))
  rawKern <- t(matrix(dot$K, nc, nr))
  if (is.null(theta)) theta <- 1/max(rawKern)
  K <- exp(-rawKern * theta)
  class(K) <- "kern"
  object <- list(K = K, theta = theta)
 #class(object) <- "gaussKern"
  return(object)
}
