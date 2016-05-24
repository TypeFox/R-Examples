NULL


#' 
#' Extracts the elements of a data frame corresponding to a period between \code{year_min} and \code{year_max}  for the stations listed in \code{station}
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'
#' @param data a dataframe containing daily data. 
#' @param year_min start year
#' @param year_max end year
#' @param station character vector of the IDs of the station where the data are required
#' 

#' @export 
#'   
#'
#' @return   a matrix containing the requested daily data where each day corresponds to a row and each station corresponds to a column
#'      
#' 
#' 
#' 
#' 

#' @note The input data frame \code{data} must have the following fields: \code{year,month,day,variables_ID1,variables_ID2,...} 
#' where the fields \code{,variables_ID1,variables_ID2,...} contain the daily variables referred to the respective stations and the field names are replaced with the respective station ID. 



extractyears <-
function(data,year_min=1961,year_max=1990,station=c("T0001","T0014","T0129"))  {
	
	
	
	station <- station[station %in% names(data)] 
	
	nstation=length(station)	
	iyear <- ((data$year<=year_max) & (data$year>=year_min))
	out <- data[iyear,station]
# See correction on getMOnthlyMean EC 2013-04-24
#	out <- subset(data,year<=year_max & year>=year_min,select=station[1])
	
#	if (nstation>1) {
#		
#		
#		for (c in 2:nstation) {
#			temp <- subset(data,year<=year_max & year>=year_min,select=station[c])
#			
#			out <- cbind(out,temp)
			
#		}
#	}
	
	return(as.matrix(out))
	
}

