#' Take a variable bounded above/below/both and return an unbounded (normalized) variable.
#'
#' This transforms bounded variables so that they are not bounded.
#' First variables are coerced away from the boundaries. by a distance of \code{tol}.
#' The natural log is used for variables bounded either above or below but not both.
#' The inverse of the standard normal cumulative distribution function 
#'   (the quantile function) is used for variables bounded above and below.
#'
#' @param x A vector, matrix, array, or dataframe with value to be 
#'   coerced into a range or set.
#' @param constraints A list of constraints.  See the examples below 
#'   for formatting details.
#' @param tol Variables will be forced to be at least this far away 
#'   from the boundaries.
#' @param trim If TRUE values in x < lower and values in x > upper 
#'   will be set to lower and upper, respectively, before normalizing.
#' @return An object of the same class as \code{x} with the values 
#'   transformed so that they spread out over any part of the real 
#'   line.
#'
#' A variable \code{x} that is bounded below by \code{lower} is
#'   transformed to \code{log(x - lower)}.
#'
#' A variable \code{x} that is bounded above by \code{upper} is
#'   transformed to \code{log(upper - x)}.
#'
#' A variable \code{x} that is bounded below by \code{lower} and
#'   above by \code{upper} is transformed to 
#'   \code{qnorm((x-lower)/(upper - lower))}.
#' @export
#' @examples
#'   constraints=list(lower=5)           # lower bound when constrining to an interval
#'   constraints=list(upper=10)          # upper bound when constraining to an interval
#'   constraints=list(lower=5, upper=10) # both lower and upper bounds
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
NormalizeBoundedVariable <- 
function(x, 
  constraints,
  tol=stats::pnorm(-5),
  trim=TRUE
) {
  if( is.null(constraints$lower) ) constraints$lower <- -Inf
  if( is.null(constraints$upper) ) constraints$upper <- Inf
  if( constraints$upper < constraints$lower ) stop("'upper' must be greater than 'lower.'")
  if( trim ) {
    x <- pmax(x, constraints$lower)
    x <- pmin(x, constraints$upper)
  } else {
    if( min(x) < constraints$lower ) stop("All values in x must be greater than or equal to the lower bound.")
    if( max(x) < constraints$upper ) stop("All values in x must be less than or equal to the upper bound.")
  }
  
  # force values away from boundaries
  if( is.finite(constraints$lower) ) x <- pmax(constraints$lower + tol, x)
  if( is.finite(constraints$upper) ) x <- pmin(constraints$upper - tol, x)
  
  if( is.infinite(constraints$lower) && is.infinite(constraints$upper) ) {
    # not bounded; degenerate case
    return( x )
  } else if( is.infinite(constraints$lower) ) {
    # only bounded above
    return( log(constraints$upper - x) )
  } else if( is.infinite(constraints$upper) ) {
    # only bounded below
    return( log(x - constraints$lower) )
  } else {
    # bounded above and below
    return( stats::qnorm((x-constraints$lower)/(constraints$upper-constraints$lower)) )
  }
}
