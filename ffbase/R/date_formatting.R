#' Date Conversion Functions for \code{ff} vector
#' 
#' Date Conversion Functions for \code{ff} vector.
#'
#' @export
#' @export as.Date.ff_vector
#' @method as.Date ff_vector
#' @param x an object of class \code{ff_vector}
#' @param ... other parameters passed on to \code{\link{as.Date}}
#' @param inplace passed on to \code{\link{chunkify}}
#' @return An \code{ff_vector} of length(x) containing the result of as.Date applied to the elements in chunks
as.Date.ff_vector <- chunkify(fun = as.Date)

#' Date Conversion Functions for \code{ff} vector
#' 
#' Date Conversion Functions for \code{ff} vector.
#'
#' @export
#' @export format.ff_vector
#' @method format ff_vector
#' @param x an object of class \code{ff_vector}
#' @param ... other parameters passed on to \code{\link{format}}
#' @param inplace passed on to \code{\link{chunkify}}
#' @return An \code{ff_vector} of length(x) containing the result of format applied to the elements in chunks
format.ff_vector <- chunkify(fun = format)



