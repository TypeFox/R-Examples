#' @export which.nearest
#' 
#' @title Which Nearest
#' @description Find values of one vector that are nearest to values in another vector.
#' 
#' @param x vector of values to be compared against.
#' @param y vector of values to examine relative to \code{x}. May be of length 1.
#' 
#'  @return For each value in \code{y}, returns index of value of \code{x} which is 
#'    nearest to \code{y} in absolute value. In the case of ties, the function returns the
#'    first index of \code{x}. If nearest value is \code{min(x)} or \code{max(x)}, a warning is issued. 
#'    \code{NA}s and \code{NaN}s in \code{x} are ignored; \code{NA}s and \code{NaN}s in \code{y} are returned.
#'    
#' @author Tim Gerrodette \email{tim.gerrodette@@noaa.gov}
#' 
#' @examples
#' x <- sort(sample(1:100, 20))
#' y <- sort(sample(min(x):max(x), 5))
#' i <- which.nearest(x, y)
#' x
#' y
#' x[i]

which.nearest <- function(x, y) {
  sapply(y, function(i) {
    if (is.na(i) | is.nan(i)) return(i)
    if (i < min(x, na.rm = TRUE)) warning("y = ", i, " is less than minimum of x", call. = FALSE)
    if (i > max(x, na.rm = TRUE)) warning("y = ", i, " is greater than maximum of x", call. = FALSE)
    return(which.min(abs(i - x)))
  })
}