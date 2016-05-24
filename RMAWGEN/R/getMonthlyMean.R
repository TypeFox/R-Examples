
NULL

#' 
#' Calculates the monthly means of a data frame corresponding to a period between \code{year_min} and \code{year_max}  for stations listed in \code{station}
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'
#' @param data a dataframe containing daily data. 
#' @param year_min start year
#' @param year_max end year
#' @param station character vector of the IDs of the station where the data are requested
#' @param no_date logical value if \code{TRUE} the function \code{extractmonths} is used. 
#' Default is \code{FALSE}. It is recommended if \code{data} does not contain columns for the dates. 
#' @param yearly  logical value. If \code{TRUE} the monthly mean values are calculated for each year from \code{year_min} to \code{year_max} separately. Default is \code{FALSE}.
#' @param origin date corresponding to the first row
#'     
#' @export  
#'
#' @return   a matrix containing the requested monthly means where each month corresponds to a row and each station corresponds to a column or a list of such matrices in case the monthly mean values are calculated separately for each year (if \code{yearly} is \code{TRUE})
#'      
#' @import date
#' 
#' @seealso \code{\link{extractyears}}
#' 

#' @note The input data frame \code{data} must have the following fields: \code{year,month,day,variables_ID1,variables_ID2,...} 
#' where the fields \code{,variables_ID1,variables_ID2,...} contain the daily variables referred to the respective stations and the field names are replaced with the respective station ID. In case \code{yearly} is \code{TRUE} the returned output is a list of matrices whose names are the corresponding year. 

# temp station=c("T0001","T0014","T0129")

getMonthlyMean <-
function(data,year_min=1961,year_max=1990,station=names(data),no_date=FALSE,origin="1961-1-1",yearly=FALSE)  {
	
	

	if (is.null(station)) {
	
		nstation=ncol(data) 
	
	} else { 		
	
		nstation=length(station)
	
	}



	if (no_date) {
		
		dates <- findDate(1:nrow(data),origin=origin,data.frame=TRUE,decimal=FALSE,character=FALSE)
		
		data <- cbind(dates,data) 
		datalocked <- data
		
	}	

	if (yearly) {
		out <- list()
		
		for (year in year_min:year_max) {
			out[[as.character(year)]] <- getMonthlyMean(data,year_min=year,year_max=year,station=station,no_date=FALSE,origin=NULL,yearly=FALSE)	
			
		}
		
		
		
		
	} else {
		
		nmonth=12
	
	
		out <- as.data.frame(array(NA,c(nmonth,nstation)))
		names(out) <- station
		
		for (r in 1:nmonth) {
			im <- (data$year<=year_max) & (data$year>=year_min) & (data$month==r)
			for (c in 1:nstation) {
				var <- data[im,station[c]]
			
				
	# The above line repleces according to 'subset' help page warnings: var <- subset(data,year<=year_max & year>=year_min & month==r,select=station[c])
		
				out[r,c] <- mean(var[!is.na(var)])
			}
		}
		
	}

	

	return(out)
	
}

