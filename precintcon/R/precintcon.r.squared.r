#' @noRd
#' @name precintcon.r.squared
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.r.squared 
#' @title Coefficient of determination
#' @description Calculates the coefficient of determination for the concentration
#' index model. 
#' @usage precintcon.r.squared(object, a, b) 
#' @param object is a daily precipitation serie.
#' @param a,b are constants defined in the concentration index.
#' @return The coefficient of determination. 
#' @seealso 
#' \code{\link{ci}}
#' \code{\link{read.data}} 
#' @keywords precipitation R squared 
precintcon.r.squared <- function(object, a, b) {
	if (is.element("precintcon.fd", class(object))) {
		y.tilde <- a * object[,8] * exp(b * object[,8])
		return( 1 - (sum((object[,9]-y.tilde)^2)/sum((object[,9]-mean(object[,9]))^2)) )

	} else
		stop("Invalid data. Please, check your input object.")
}