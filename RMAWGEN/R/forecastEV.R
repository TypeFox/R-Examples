
NULL


#' 
#' Forecasts the expected value of a VAR realization given the prievious one 
#' 
#' @param var A VAR model represented by a \code{varest} object as returned by \code{\link{getVARmodel}} or \code{\link{VAR}}
#' @param xprev previous status of the random variable
#' @param exogen vector containing the values of the "exogen" variables (predictor) for the generation
#' 
#'  @seealso \code{\link{forecastResidual}}
#'  @author  Emanuele Cordano, Emanuele Eccel
#'    
#'  @export
#'
#'      
#' 
#' 
#' @return  a vector of values





forecastEV <-
function(var,xprev=NULL,exogen=NULL) {

	
#	class(var) <- "varest"

	if (is.null(xprev)) xprev <- rnorm(var$p*var$K)

	if (!is.null(exogen)) xprev <- c(xprev,exogen)

	cvar <- coef(var)
	out <- array(NA,length(cvar))
	for (i in 1:length(out)) {
	
		out[i] <- as.double(cvar[[i]][,"Estimate"] %*% xprev)
		
		
	}
	
	
	return(out)
	
}

