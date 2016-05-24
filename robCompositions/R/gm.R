#' geometric mean
#' 
#' This function calculates the geometric mean.
#' 
#' Calculates the geometric mean of all positive entries of a vector.
#' 
#' @param x A numeric vector.
#' @return The geometric mean.
#' @author Matthias Templ
#' @keywords math
#' @export
#' @examples
#' 
#' gm(runif(100))
#' 
gm <- function (x) {
	if(!is.numeric(x)) stop("x has to be a vector of class numeric")
	if (any(na.omit(x == 0)))
		0
	else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}
