##' Enforce array and convert to vector if appropriate
##'
##' \code{ensuredim} turns vectors into 1d-arrays, and leaves arrays unchanged. \code{drop1d} is the
##' inverse: it converts 1d arrays into vectors.
##' @param x vector (or array)
##' @return \code{esuredim} array of at least one dimension
##' @author Claudia Beleites
##' @export 
ensuredim <- function (x){
  if (is.null (dim (x)))
    x <- structure (x, .Dim = length (x),
                    .Dimnames =  lon (list (names (x))),
                    .Names = NULL)
  x
}

.test (ensuredim) <- function (){
  checkEquals (ensuredim (v), structure(1:3, .Dim = 3L, .Dimnames = list(c("a", "b", "c"))))
  checkEquals (ensuredim (as.numeric (v)), structure(1:3, .Dim = 3L))
  checkEquals (ensuredim (m), m)
  checkEquals (ensuredim (a), a)
}

##' @rdname ensuredim
##' @param drop if \code{TRUE}, 1d arrays are converted into vectors
##' @return \code{drop1d} vector, if \code{x} had only 1 dimension
##' @author Claudia Beleites
##' @export 
drop1d <- function (x, drop = TRUE){
  if (drop && ndim (x) == 1)
    x <- structure (x, .Dim = NULL, .Dimnames = NULL, .Names = dimnames (x)[[1]])

  x
}

.test (drop1d) <- function () {
  checkEquals (drop1d (a), a)
  checkEquals (drop1d (m), m)
  checkEquals (drop1d (v), v)
  checkEquals (drop1d (ensuredim (v)), v)
  checkEquals (drop1d (ensuredim (v), drop = FALSE), ensuredim (v))
}

