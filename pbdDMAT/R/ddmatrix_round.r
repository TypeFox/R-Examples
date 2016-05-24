#' Rounding of Numbers
#' 
#' Extensions of R rounding functions for distributed matrices.
#' 
#' Rounding to a negative number of digits means rounding to a power of ten, so
#' for example \code{round(x, digits = -2)} rounds to the nearest hundred.
#' 
#' @param x 
#' numeric distributed matrix
#' @param digits 
#' integer indicating the number of decimal places
#' (\code{round()}) or significant digits (\code{signif()}) to be used.
#' Negative values are allowed (see 'Details').
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
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- ceiling(x/3)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name rounding
#' @rdname rounding
NULL

#' @rdname rounding
#' @export
setMethod("round", signature(x="ddmatrix"),
  function(x, digits=0)
  {
    x@Data <- round(x@Data, digits=digits)
    return(x)
  }
)

#' @rdname rounding
#' @export
setMethod("ceiling", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- ceiling(x@Data)
    return(x)
  }
)

#' @rdname rounding
#' @export
setMethod("floor", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- floor(x@Data)
    return(x)
  }
)

