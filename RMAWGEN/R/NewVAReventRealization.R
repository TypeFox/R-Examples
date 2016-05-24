
NULL


 
#'   
#' Generates a new realization of a VAR model 
#'  
#'  
#'
#' @param var A VAR model represented by a \code{varest} object as returned by \code{\link{getVARmodel}} or \code{\link{VAR}}
#' @param xprev previous status of the random variable
#' @param noise uncorrelated or white noise (residual). Default is \code{rnorm(length(xprev))} (or \code{rnorm(ncol(B)})
#' @param exogen vector containing the values of the "exogen" variables (predictor) for the generation
#' @param B matrix of coefficients for the vectorial white-noise component
#' 
#'       
#' 
#' @export 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' @seealso \code{\link{forecastEV}},\code{\link{forecastResidual}}
#' 
#' 
#' 
#' @return  a vector of values




NewVAReventRealization <-
function(var,xprev,noise,exogen=NULL,B=NULL) {
	
	
	
	out <- NULL 
	
	
	out <- forecastEV(var=var,xprev=xprev,exogen=exogen)+forecastResidual(var=var,xprev=noise,B=B)
	

	
	
	return(out)
	
}

