# Created by Matthew A. Birk
# Calculates index of closest value
# Last updated: Jan 2016

#' Where is the closest?
#'
#' Returns index of the closest value to \code{x}.
#'
#' @param vec a numeric vector.
#' @param x numeric. The value for which the closest match should be returned.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @seealso \code{\link{which.min}}, \code{\link{which.max}}
#'
#' @examples
#' which.closest(10:1, 3.3)
#' 
#' @encoding UTF-8
#' @export

which.closest = function(vec, x){
	which.min(abs(vec - x))
}