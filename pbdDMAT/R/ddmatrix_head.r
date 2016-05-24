#' Head and Tail of a Distributed Matrix
#' 
#' The functions rely on reblocking across different BLACS contexts.  If
#' \code{i} is not empty, then the input distributed matrix will be
#' redistributed along context 1, where extracting/deleting rows does not
#' destroy block-cyclicality. Likewise, if \code{j} is not empty, then the
#' input distributed matrix will be redistributed along context 2. When
#' extraction is complete, the matrix will be redistributed across its input
#' context.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param n 
#' a single integer. If positive, size for the resulting object:
#' number of elements for a vector (including lists), rows for a matrix or data
#' frame or lines for a function. If negative, all but the \code{n} last/first
#' number of elements of \code{x}.
#' @param ... 
#' additional arguments.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:9, 3, bldim=2)
#' 
#' y <- head(y, 2)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction
#' @name headsortails
#' @rdname headsortails
NULL


headsortails <- function(x, n, index)
{
  n <- as.integer(n)
  dim <- as.integer(dim(x)[1])
  
  if (n == 0 || (n < 0 && -n>dim) )
    return(x[0, ])
  else if (n < 0)
  {
    n <- dim+n
    return(x[index, ])
  }
  else
  {
    if (n >= dim)
      return(x)
    else
    return(x[index, ])
  }
}



#' @rdname headsortails
#' @export
head.ddmatrix <- function(x, n=6L, ...)
{
  index <- 1L:n
  headsortails(x=x, n=n, index=index)
}



#' @rdname headsortails
#' @export
tail.ddmatrix <- function(x, n=6L, ...)
{
  index <- (dim-n+1L):n
  headsortails(x=x, n=n, index=index)
}

