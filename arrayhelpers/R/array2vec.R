##' Converting array and vector Indices
##' Calculate the vector index from array indices, and vice versa.
##' 
##' \code{array}s are \code{numeric}s with a \code{dim} attribute and are
##' stored with the first index moving fastest (i.e. by column). They can be
##' indexed both ways.
##' 
##' @aliases array2vec vec2array
##' @param iarr vector with the indices into the array dimensions
##' @param dim vector with the array dimensions, as returned by \code{dim (x)}
##' @return \code{array2vec} returns a scalar, \code{vec2array} a
##'   \code{matrix}.
##' @author C. Beleites
##' @seealso see \code{\link[base]{Extract}} on the difference of indexing an
##'   \code{array} with a vector or a \code{matrix}.
##' @keywords array
##' @export
##' @examples
##' 
##' arr <- array (rnorm (24), dim = 2 : 4)
##' arr
##' 
##' v <- matrix(c(2, 2, 2), nrow = 1)
##' i <- array2vec (v, dim = dim (arr))
##' i
##' arr[v]
##' arr[i]
##' 
##' arr[c(2, 2, 2)] ## indexing with a vector
##' arr[2]
##' 
array2vec <- function (iarr, dim){
  if (!is.matrix (iarr))
    dim (iarr) <- c(1, length (iarr))

  if (ncol (iarr) != length (dim))
    stop ("Number of columns in iarr and number of dimensions differ.")

  if (any (sweep (iarr, 2, dim) > 0))
    stop ("array index > dim")

  pdim <- c(1, cumprod (dim [- length (dim)]))
  iarr <- iarr - 1

  colSums(apply (iarr, 1, "*", pdim)) + 1
}


##' @param ivec scalar with the index into the vector
##' @rdname array2vec
##' @export
##' @examples
##'  
##' i <- 14
##' v <- vec2array (i, dim = dim (arr))
##' v
##' arr [v]
##' arr [i]
##' 
vec2array <- function (ivec, dim) {
  ndim <- length (dim)
  pdim <- c(1, cumprod (dim))

  iarr <- matrix(NA, nrow = length(ivec), ncol = ndim) # matrix for the array indices
  colnames (iarr) <- letters[8 + seq_len (ndim)]       # i, j, k, ...

  ivec <- (ivec - 1)
  for (j in seq_len (ndim))
    iarr [, j] <- (ivec %% pdim [j + 1]) / pdim [j]

  1 + floor(iarr)
}
