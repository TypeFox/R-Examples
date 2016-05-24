
NULL


#' 
#' Forecasts the residual value of a VAR realization given the white noise covariance matrix
#' 
#' 
#' @param var A VAR model represented by a \code{varest} object as returned by \code{\link{getVARmodel}} or \code{\link{VAR}}
#' @param xprev previous status of the random variable, in this case the "current instant"white-noise". Default is \code{NULL} and then randomly generated.
#' @param B matrix of coefficients for the vectorial white-noise component
#' 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'    
#' 
#' @seealso \code{\link{forecastEV}},\code{\link{NewVAReventRealization}}
#' 
#'        
#' 
#' @export
#' 
#' @return  a vector of values




forecastResidual <-
function(var,xprev=NULL,B=NULL) {

	
	if (is.null(B)) B <- t(chol(summary(var)$covres))
	
	if (is.null(xprev)) xprev <- rnorm(ncol(B))
	
	
	out <- as.vector(B %*% xprev)
	
	return(out)
	
}

