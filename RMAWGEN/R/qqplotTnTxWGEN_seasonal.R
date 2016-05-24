
NULL

#'
#'  Makes four seasonal qqplots  (winter, spring, summer and autumn) of measured and simulated data for several stations.   
#' 
#'  
#' 
#' @param measured  matrix containing measured data (each station corresponds to a column) 
#' @param simulated matrix containing respective generated data (each station corresponds to a column) 
#' @param xlab,ylab see \code{\link{plot.default}},\code{\link{qqplotWGEN}}
#' @param titles titles that will be added 
#' @param station character vector containing IDs of analyzed station. If \code{NULL} (default) all station (columns of \code{simulated} and \code{measured}) are considered
#' @param directorypdf name of the directory (path included) where to seva the outputs 
#' @param origin first day of data, see \code{\link{extractmonths}} for format and other information 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#' @export
#' 
#' 
#' @note  Uses \code{\link{qqplotTnTxWGEN}} for each seasons of collected data and saves the output on  pdf files. See the R code for further details.
#' @seealso \code{\link{qqplotTnTxWGEN}},\code{\link{extractmonths}}
#' 
#' 
#' 
#'        
#' @return  0 in case of success 



qqplotTnTxWGEN_seasonal <-
function (measured,simulated,origin="1961-1-1",xlab="simulated[degC]",ylab="measured[degC]",titles=c("Q-Qplot_An._Tx","Q-Qplot_An._Tn"),directorypdf,station=NULL) {
	
	
	
	winter <- c("Dec","Jan","Feb")
	spring <- c("Mar","Apr","May")
	summer <- c("Jun","Jul","Aug")
	autumn <- c("Sep","Oct","Nov")
	
	mes_winter <- extractmonths(data=measured,origin=origin,when=winter) 
	mes_spring <- extractmonths(data=measured,origin=origin,when=spring)
	mes_summer <- extractmonths(data=measured,origin=origin,when=summer)
	mes_autumn <- extractmonths(data=measured,origin=origin,when=autumn)
	
	sim_winter <- extractmonths(data=simulated,origin=origin,when=winter) 
	sim_spring <- extractmonths(data=simulated,origin=origin,when=spring)
	sim_summer <- extractmonths(data=simulated,origin=origin,when=summer)
	sim_autumn <- extractmonths(data=simulated,origin=origin,when=autumn)
	
	filename_winter <- paste(directorypdf,"winter.pdf",sep="/")
	filename_spring <- paste(directorypdf,"spring.pdf",sep="/")
	filename_summer <- paste(directorypdf,"summer.pdf",sep="/")
	filename_autumn <- paste(directorypdf,"autumn.pdf",sep="/")
	
	pdf(filename_winter)
	qqplotTnTxWGEN(measured=mes_winter,simulated=sim_winter,xlab=xlab,ylab=ylab,titles=titles,station=station)
	dev.off()
	
	pdf(filename_spring)
	qqplotTnTxWGEN(measured=mes_spring,simulated=sim_spring,xlab=xlab,ylab=ylab,titles=titles,station=station)
	dev.off()
	
	pdf(filename_summer)
	qqplotTnTxWGEN(measured=mes_summer,simulated=sim_summer,xlab=xlab,ylab=ylab,titles=titles,station=station)
	dev.off()
	
	pdf(filename_autumn)
	qqplotTnTxWGEN(measured=mes_autumn,simulated=sim_autumn,xlab=xlab,ylab=ylab,titles=titles,station=station)
	dev.off()
	
	
	
	
	return(0)
}

