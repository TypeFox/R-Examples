##' number of dimensions
##'
##' @param a vector, matrix, or array
##' @param ... indexing instructions. The names of the arguments specify the dimension 
##'    (i = 1st, j = 2nd, ...). The indexing expressions are the same as for \code{\link[base]{[}}
##' @return integer: length of dim attribute
##' @author Claudia Beleites
##' @export 
ndim <- function (a){
  length (dim (a))
}

.test (ndim) <- function (){
  checkEquals (ndim (v), 0L)
  checkEquals (ndim (ensuredim (v)), 1L)
  checkEquals (ndim (m), 2L)
  checkEquals (ndim (a), 3L)
}
