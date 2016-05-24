##' Transpose arrays
##'
##' This function provides transposing of arrays or vectors as swapping their first two dimensions.
##' \code{t (array)} can be enabled via \code{\link[methods]{setMethod}}, see the example.
##' @param x an array 
##' @return the array with the first two dimensions swapped.
##' @author Claudia Beleites
##' @seealso \code{\link[base]{t}}
##' @export 
##' @examples
##' a <- array (1 : 24, 4:2)
##' a
##' ta (a)
##' 
##' setMethod ("t", "array", ta)
##' t (a)
##' removeMethod ("t", "array")
##'
ta <- function (x){
  if (! (is.vector (x) || is.matrix (x) || is.array (x)))
      stop ("x must be array, matrix, or vector.")
      
  a <- makeNd (x, 2)                 # ensure at least 2 dimensions

  d <- seq_along (dim (x))
  d [1 : 2] <- 2 : 1                    # swap first 2 dimensins
  
  aperm (x, d)
}

.test (ta) <- function () {
  checkEqualsNumeric (a, ta (ta (a)))
  checkEquals (a, ta (ta (a)))
  checkIdentical (a, ta (ta (a)))


  checkIdentical (ta (m), t (m))

#  instd <- length (findMethod ("t", "array")) > 0

#  if (instd) removeMethod ("t", "array")
#  checkException (t (a))
#  setMethod ("t", "array", ta)
#  checkIdentical (t (t (a)), a)
  
#  if (!instd) removeMethod ("t", "array")
}

