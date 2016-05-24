
NULL

 
#'
#'  Makes a qqplot of measured and simulated data for several stations. 
#' 
#'  
#' @param measured  matrix containing measured data (each station corresponds to a column) 
#' @param simulated matrix containing respective generated data (each station corresponds to a column) 
#' @param xlab,ylab see \code{\link{plot.default}},\code{\link{qqplotWGEN}}
#' @param titles titles that will be added to \code{main} argument of \code{\link{plot.default}}
#' @param station character vector containing IDs of analyzed station. If \code{NULL} (default) all station (columns of \code{simulated} and \code{measured}) are considered
#' @param diff,quantile see \code{\link{qqplotWGEN}}
#' 
#' @note It uses \code{\link{qqplotWGEN}} and makes a figure for each pair of columns from \code{measured} and \code{simulated}. See the R code for further details.
#' 
#' @export
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#'     
#' @return  0 in case of success 
#' 
#' 








qqplotTnTxWGEN <-
function (measured,simulated,xlab="simulated[degC]",ylab="measured[degC]",titles=c("Q-Qplot_An._Tx","Q-Qplot_An._Tn"),station=NULL,diff=FALSE,quantile=0) {
	

	
	
	nrows=2
	ncols=ncol(measured)/2
#	str(measured)
	name <- names(simulated)
	if (is.null(station)) {
		ncols=ncol(measured)/2
		istation=1:ncols
	} else {
		
		istation=which(name %in% station)
#		str(istation)
		ncols=length(istation)/2
	}
	
	
	
	#print(nrows)
	#print(ncols)
	
	par(mfrow=c(nrows,ncols))
	
	
	ir=0
	xlim <- c(min(simulated,na.rm=TRUE),max(simulated,na.rm=TRUE))
	ylim  <- c(min(measured,na.rm=TRUE),max(measured,na.rm=TRUE))
	for (r in 1:nrows) {
		for (c in istation[1:ncols]) { 	
		
			main=paste(titles[r],name[c],sep=" ")
			qqplotWGEN(val=cbind(simulated[,c+ir],measured[,c+ir]),ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,main=main,diff=diff,quantile=quantile)
		}
		ir=ir+ncols
	}
	
#	dev.off()
	return(0)	
}

