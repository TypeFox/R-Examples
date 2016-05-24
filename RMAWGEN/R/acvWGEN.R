

NULL


#'  
#' Plots the auto- and cross- covariance functions between measured and simulated data for several stations
#' 
#' 
#' 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'    
#'  
#'
#' @param  measured  matrix containing measured time series
#' @param  simulated matrix containing simulated time series 
#' @param  titles  title suffixes for the simulated and measured data respectively c("Sim.","Mes.")
#' @param  station string vector containing the IDs of the meteorological stations where the autocovariance is calculated. 
#' If it is \code{NULL} (default) all stations (corresponding to the columns of "simulated" and "measured") are applied
#' 
#'  
#' @export 
#' @note It uses \code{\link{acf}} function
#' @return  0 in case of success 

#xlab="simulated[° C]",ylab="measured[° C]"
#


acvWGEN <-
function (measured,simulated,titles=c("Sim.","Mes."),station=NULL) {

	
	val <- cbind(simulated,measured)
	if (!is.null(station)) {
		
		val <- val[,names(val) %in% station]
	}
	
	
	nrows=2
	ncols=ncol(val)/2
	
	
	
	
	name <- names(val)
	
	par(mfrow=c(nrows,ncols))
	
	
	ir=0
	
	
	for (r in 1:nrows) {
		
		for (c in 1:ncols) { 	
			
			main=paste(titles[r],name[c],sep=" ")
			acf(val[,c+ir],main=main)	
		}
		ir=ir+ncols
	}
	
	
	return(0)	
}

