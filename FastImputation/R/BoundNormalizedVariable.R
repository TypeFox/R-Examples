#' Take a normalized variable and transform it back to a bounded variable.
#'
#' This takes variables on the real line and constrains them to be on
#' a half-line (constrained above or below) or a segment (constrained both
#' above and below). This is approximately the inverse of 
#' \code{NormalizeBoundedVariable}; this does not completely reverse the
#' effect of \code{NormalizeBoundedVariable} because \code{NormalizeBoundedVariable}
#' first forces values away from the bounds, and this information is lost.
#'
#' @param x A vector, matrix, array, or dataframe with value to be coerced into a range or set.
#' @param constraints A list of constraints.  See the examples below for formatting details.
#' @return An object of the same class as x with the values transformed into the desired half-line or segment.
#' @export
#' @examples
#'   constraints=list(lower=5)           # lower bound when constrining to an interval
#'   constraints=list(upper=10)          # upper bound when constraining to an interval
#'   constraints=list(lower=5, upper=10) # both lower and upper bounds
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
BoundNormalizedVariable <- 
function(x, 
  constraints
) {
  if( is.data.frame(x) ) {
    was.data.frame <- TRUE
    df.names <- names(x)
    x <- as.matrix(x)
  } else {
    was.data.frame <- FALSE
  }
  if( !is.numeric(x) ) stop("x must be numeric")
  if( !(is.vector(x) || is.matrix(x) || is.array(x)) ) stop("x must be a vector, matrix, or array")
  if( is.null(constraints$lower) ) constraints$lower <- -Inf
  if( is.null(constraints$upper) ) constraints$upper <- Inf
  if( constraints$upper < constraints$lower ) stop("'upper' must be greater than 'lower.'")
  
  if( is.infinite(constraints$lower) && is.infinite(constraints$upper) ) {
    # not bounded; degenerate case
    out <- x
  } else if( is.infinite(constraints$lower) ) {
    # only bounded above
    out <- constraints$upper - exp(x)
  } else if( is.infinite(constraints$upper) ) {
    # only bounded below
    out <- exp(x) + constraints$lower
  } else {
    # bounded above and below
    out <- stats::pnorm(x) * (constraints$upper-constraints$lower) + constraints$lower
  }
  if( was.data.frame ) {
    out <- as.data.frame(out)
    names(out) <- df.names
  }
  return( out )
}
