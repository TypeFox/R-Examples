
NULL

#' 
#' Calculates the daily means of a range of days around each date of a data frame corresponding to a period between \code{year_min} and \code{year_max}  for stations listed in \code{station}
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'
#' @param data a data frame containing daily data. 
#' @param year_min start year
#' @param year_max end year
#' @param station character vector of the IDs of the station where the data are requested

#' @param origin origin date of time-series
#' @param lag lag (number of days) on which daily mean is calculated. The mean is calculated considereing \code{lag} days before and after each day.  
#'
#' @return  a matrix containing the requested daily mean data where each day corresponds to a row and each station corresponds to a column
#'      
#' 
#' 
#' @seealso \code{\link{extractyears}}
#' 
#' 
#' @export
#' 
#' @note The input data frame \code{data} must have the following fields: \code{year,month,day,variables_ID1,variables_ID2,...} 
#' where the fields \code{,variables_ID1,variables_ID2,...} contain the daily variables referred to the respective stations and the field names are replaced with the respective station ID. 

# @param no_date logical value. If \code{TRUE} the function \code{extractmonths} is used, it is recommended if \code{data} does not contain columns for the dates. 
# Default is \code{FALSE}


getDailyMean <-
function(data,year_min=1961,year_max=1990,station=c("T0001","T0010"),origin="1961-1-1",lag=5) {


	
	
	dates <- findDate(1:nrow(data),origin=origin,data.frame=FALSE,decimal=FALSE,character=FALSE)
	iday <-  as.integer(format(dates,"%j"))
	year_v <-  as.integer(format(dates,"%Y"))
	iday2 <- iday[(year_v>=year_min) & (year_v<=year_max)] # cutoff iday between year_min and year_max 
	
	coln <- length(station)
	rown <- nrow(data[year_v>=year_min & year_v<=year_max,])
	
	out <- as.data.frame(array(NA,c(rown,coln)))
	ndoy<- 365
	for (c in 1:coln) {
		for (r in 1:rown) {
			
			
			out[r,c] <- mean(data[(year_v>=year_min) & (year_v<=year_max) & 
									(((iday>=iday2[r]-lag) & (iday<=iday2[r]+lag) ) | 
									((iday+ndoy>=iday2[r]-lag) & (iday+ndoy<=iday2[r]+lag) ) |
									((iday-ndoy>=iday2[r]-lag) & (iday-ndoy<=iday2[r]+lag) )),station[c]],na.rm=TRUE)
			 
		}
	}
	
 

	
	names(out) <- station

	return(out)
}

