# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL
#'
#' It calculates the number of wet days for each month and each year
#' 
#' @param data data frame R object containing daily precipitation time series for several gauges (one gauge time series per column). 
#' @param valmin threshold precipitation value [mm] for wet/dry day indicator.
#' @param origin character string \code{"yyyy-mm-dd"} indicated the date of the first row of \code{"data"}. 
#' @param station character string indicating the stations. Default is \code{names(data)}
#' @export
#'
#' 
#' @return Function returns a list of data frames containing the spell length expressed in days
#' 
#' @examples  
#' 
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' 
#' period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' prec_mes <- PRECIPITATION[period,station]  
#' 
#' ## removing nonworking stations (e.g. time series with NA)
#' accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#' for (it in names(prec_mes)) {
#' 		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#' }
#'
#' prec_mes <- prec_mes[,accepted]
#' ## the dateset is reduced!!! 
#' prec_mes <- prec_mes[,1:3]
#' 
#' origin <- paste(year_min,1,1,sep="-")
#' 
#' nwetdays <- nwetdays(prec_mes,origin)
#' 


nwetdays <- function(data,valmin=0.5,origin="1961-1-1",station=names(data)) {
	
	station <- station
	data <- adddate(data,origin=origin)
	
	months <- unique(data$month)
	years <- unique(data$year)
	out <- as.data.frame(array(0,c(length(years)*length(months),length(station)+2)))
	
	names(out) <- c("month","year",station) 
	str(station)
	str(out)
	out$year <- rep(years,each=length(months))
	out$month <- rep(months,times=length(years))
	
	miyd <- sprintf("%04d-%02d",data$year,data$month)
	miyo <- sprintf("%04d-%02d",out$year,out$month)
	for (it in miyo) {
		
		temp <- data[miyd==it,station]>=valmin
		
		temp <- as.matrix(temp)
		
		row <- which(miyo==it)
		
		out[row,station] <- apply(X=temp,MARGIN=2,FUN=function(x) { length(which(x==TRUE))})
		
		
	}
	
	
	
	return(out)
	
}

	
