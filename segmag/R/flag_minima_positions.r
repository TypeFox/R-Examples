#' @rdname flag_maxima_positions
#' @export
flag_minima_positions <- function(values)
{
  if (! is.numeric(values) || length(values)==0) stop("values must be a numeric vector and be of length > 0")
  if (max(values) == min(values)) stop("There is no variation in vector values: max and min must not be equal")
  
  # Minima are the maxima of values * -1
  return( flag_maxima_positions_impl(values * -1) == 1 )
}