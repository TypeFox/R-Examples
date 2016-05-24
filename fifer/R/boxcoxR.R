#' Transform data using a boxcox transformation
#'
#' This function takes a variable as input, computes the optimal lambda
#' using a boxcox transformation, then returnes a transformed version of the variable.
#'
#' @details The \code{MASS} package has a function that computes the optimal lambda value for a particular regression equation. However, it currently (as of version 7.3-23)
#' returns a lambda vector rather than the boxcox transformed variable. This function is a wrapper for \code{\link{boxcox}} that actually returns a vector that is a 
#' transformed version of the original variable. 
#' @param x a numerical vector
#' @param minval before a transformation is performed, the variables must often be positive. This tells R what the minimum value should be. Defaults to .01. 
#' @param ... additional parameters to be used in the model fitting.
#' @export
#' @import MASS
#' @aliases boxcox.R
#' @author Dustin Fife \email{fife.dustin@@gmail.com}.
#' @examples
#' x = rnorm(100)^2
#' ### use original boxcox function
#' boxcox(x~1, plot=FALSE) ## returns a vector of lambda values and their likelihoods
#' ### use boxcoxR function
#' boxcoxR(x)
#' @references  Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#' @note This function calls the boxcox function in the MASS package. To avoid loading the package, I have branched the function directly into
#' the fifer package. 
boxcoxR = function(x, minval=.01, ...){
	variable = x
	if (min(variable, na.rm=T)<=0){
		variable = ((variable + (minval-min(variable, na.rm=T))))		#### make positive
	}
	bc = boxcox(variable~1, plot=FALSE, ...)
	lambda = bc$x[bc$y == max(bc$y)]
	if (lambda==0){
		variable = log(variable)	
	} else {
		variable = (variable^lambda-1)/lambda
	}
	return(variable)
}
