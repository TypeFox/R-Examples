#' Detect local maxima/minima of a numeric vector
#' 
#' Fast detection of local maxima and minima of a numeric vector.
#' This function takes a numeric vector as input and returns a logical vector with the
#' same length and TRUE values at local maxima/minima (depending on function).
#' If multiple succeeding values at a local maximum/minimum are equal, then only
#' the center (rounded up if necessary) of the maximum/minimum is marked with TRUE.
#' 
#' @param values numeric vector of values
#' @return logical vector with TRUE at the center of local maxima/minima
#' @examples
#' var <- c(1,2,3,3,2,1,4,5,6,7,5,4,3)
#' 
#' ## Using the Maxima functions
#' flag_maxima_positions(var)
#'
#' # values of maxima
#' var[flag_maxima_positions(var)]
#'
#' # indices of maxima
#' which(flag_maxima_positions(var))
#' 
#' 
#' ## Using the Minima functions
#' flag_minima_positions(var)
#'
#' # values of maxima
#' var[flag_minima_positions(var)]
#'
#' # indices of maxima
#' which(flag_minima_positions(var))
#' @export
flag_maxima_positions <- function(values)
{
  if (! is.numeric(values) || length(values)==0) stop("values must be a numeric vector and be of length > 0")
  if (max(values) == min(values)) stop("There is no variation in vector values: max and min must not be equal")
  
  return( flag_maxima_positions_impl(values) == 1 )
}