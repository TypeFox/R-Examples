
NULL

 
#' 
#' Makes a qqplot of measured and simulated data for several stations. 
#' 
#'  
#' @param measured  matrix containing measured data (each station corresponds to a column) 
#' @param simulated matrix containing respective generated data (each station corresponds to a column) 
#' @param xlab,ylab  see \code{\link{plot.default}},\code{\link{qqplotWGEN}}
#' @param title title 
#' @param station character vector containing IDs of analyzed stations. If \code{NULL} (default) all stations (columns of \code{simulated} and \code{measured}) are considered
#' @param diff,quantile see \code{\link{qqplotWGEN}}
#' 
#' @export
#' 
#' @note It uses \code{\link{qqplotWGEN}} and makes a figure for each pair of columns from \code{measured} and \code{simulated}. See the R code for further details.
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#'     
#' @return  0 in case of success 
#' 
#' 








qqplotprecWGEN <-
function (measured,simulated,xlab="simulated[mm]",ylab="measured[mm]",title="daily precipitation",station=NULL,diff=FALSE,quantile=0) {
	

	
	

	ncols=1 #ncol(measured)
	
	name <- names(as.data.frame(simulated))

	if (is.null(station)) {
		nrows=ncol(measured)
		istation=1:ncols
	} else {

		istation=which(name %in% station)
		nrows=length(istation)
		
	}
	
	
	
	
	
	par(mfrow=c(nrows,ncols))
	
	
	
	xlim <- c(min(simulated,na.rm=TRUE),max(simulated,na.rm=TRUE))
	ylim <- c(min(measured,na.rm=TRUE),max(measured,na.rm=TRUE))
	for (r in 1:nrows) {
		for (c in 1:ncols) { 	
		
			main=paste(title,name[r],sep=" ")
			qqplotWGEN(val=cbind(simulated[,r],measured[,r]),ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,main=main,diff=diff,quantile=quantile)
		}
	}
	

	return(0)	
}

