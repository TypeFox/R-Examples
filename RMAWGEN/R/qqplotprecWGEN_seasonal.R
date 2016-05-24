
NULL

#'
#'  Makes four seasonal qqplots  (winter, spring, summer and autumn) of measured and simulated data for several stations.   
#' 
#'  
#' 
#' @param measured  matrix containing measured data (each station corresponds to a column) 
#' @param simulated matrix containing respective generated data (each station corresponds to a column) 
#' @param xlab,ylab see \code{\link{plot.default}},\code{\link{qqplotWGEN}}
#' @param title title 
#' @param station character vector containing IDs of analyzed stations. If \code{NULL} (default) all stations (columns of \code{simulated} and \code{measured}) are considered
#' @param directorypdf name of the directory (path included) where to seva the outputs 
#' @param origin first day of data, see \code{\link{extractmonths}} for format and other information 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#' @export
#' 
#' @note  Uses \code{\link{qqplotprecWGEN}} for each season of collected data and saves the output on  pdf files. See the R code for further details.
#' @seealso \code{\link{qqplotprecWGEN}},\code{\link{extractmonths}}
#' 
#' 
#' 
#'        
#' @return  0 in case of success 



qqplotprecWGEN_seasonal <-
function (measured,simulated,origin="1961-1-1",xlab="simulated[mm]",ylab="measured[mm]",title="daily_precipitation",directorypdf,station=names(simulated)) {
	
	
	
	winter <- c("Dec","Jan","Feb")
	spring <- c("Mar","Apr","May")
	summer <- c("Jun","Jul","Aug")
	autumn <- c("Sep","Oct","Nov")
	
	mes_winter <- extractmonths(data=as.matrix(measured),origin=origin,when=winter) 
	mes_spring <- extractmonths(data=as.matrix(measured),origin=origin,when=spring)
	mes_summer <- extractmonths(data=as.matrix(measured),origin=origin,when=summer)
	mes_autumn <- extractmonths(data=as.matrix(measured),origin=origin,when=autumn)
	
	sim_winter <- extractmonths(data=as.matrix(simulated),origin=origin,when=winter) 
	sim_spring <- extractmonths(data=as.matrix(simulated),origin=origin,when=spring)
	sim_summer <- extractmonths(data=as.matrix(simulated),origin=origin,when=summer)
	sim_autumn <- extractmonths(data=as.matrix(simulated),origin=origin,when=autumn)
	
	filename_winter <- paste(directorypdf,"winter.pdf",sep="/")
	filename_spring <- paste(directorypdf,"spring.pdf",sep="/")
	filename_summer <- paste(directorypdf,"summer.pdf",sep="/")
	filename_autumn <- paste(directorypdf,"autumn.pdf",sep="/")
	
	pdf(filename_winter)
	qqplotprecWGEN(measured=mes_winter,simulated=sim_winter,xlab=xlab,ylab=ylab,title=title,station=station)
	dev.off()
	
	pdf(filename_spring)
	qqplotprecWGEN(measured=mes_spring,simulated=sim_spring,xlab=xlab,ylab=ylab,title=title,station=station)
	dev.off()
	
	pdf(filename_summer)
	qqplotprecWGEN(measured=mes_summer,simulated=sim_summer,xlab=xlab,ylab=ylab,title=title,station=station)
	dev.off()
	
	pdf(filename_autumn)
	qqplotprecWGEN(measured=mes_autumn,simulated=sim_autumn,xlab=xlab,ylab=ylab,title=title,station=station)
	dev.off()
	
	
	
	
	return(0)
}

