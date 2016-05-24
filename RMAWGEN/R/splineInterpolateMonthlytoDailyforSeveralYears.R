NULL


#' 
#' Interpolates monthly data to daily data using \code{\link{splineInterpolateMonthlytoDaily}}  for several years
#' 
#' 
#' @param val matrix containing monthly mean data for one year
#' @param start_year first year
#' @param nyear number of years since \code{start_year}
#' @param leap logical variable If \code{TRUE} (default) leap years are considered, otherwise they are not
#' @param offset integer values. Default is 2. Number of years considered beyond the extremes in order to avoid edge errors 
#' @param no_spline logical value. If \code{TRUE} no spline interpolation is calculated and the daily value corresponds to the monthly average value. Default is \code{FALSE}.
#' @param yearly logical value. If \code{TRUE} the result with men value per each month per each year. Default is \code{FALSE}.
#' 
#' 
#'  @return a matrix or data frame with interpolated daily data 
#' 
#' @export 
#' 
#' @import chron 
#' 
#' @seealso \code{\link{spline}},\code{\link{splineInterpolateMonthlytoDaily}}
#' 
#' @author Emanuele Cordano, Emanuele Eccel
#' 
#' 





splineInterpolateMonthlytoDailyforSeveralYears <-
function(val,start_year=2010,nyear=1,leap=TRUE,offset=2,no_spline=FALSE,yearly=FALSE) {
	
	
	
	
	tval <- val 
	
	
	if (yearly) {
		nyear_sim=length(val)+2*offset 
		
		tval <- val[[1]]
		for(i in 2:nyear_sim) {
			
			if (i<=offset) {
				tval <- rbind(tval,val[[1]])
			} else if (i>length(val)) {
				
				tval <- rbind(tval,val[[length(val)]])
				
			} else {
				tval <- rbind(tval,val[[i]])
			}
			
			
		}
		
		
		
	} else {
		
		nyear_sim=nyear+2*offset
		
		if (nyear_sim>1) for (i in 2:nyear_sim) tval <- rbind(tval,val)
		
		
	}
#	  
#	
#	
#	
	end_year=start_year+nyear-1
	
	count=0
	offset_day=0
	count1=0
	for (year in (start_year-offset):(end_year+offset)) {
		
		if (leap & leap.year(year)) {
		
			count=count+366
			if (year<start_year) offset_day=offset_day+366
			if (year<=end_year) count1=count1+366
		} else {

			count=count+365	
			if (year<start_year) offset_day=offset_day+365
			if (year<=end_year) count1=count1+365
		}
		
	}
#	print(offset_day) class
	nday=count
	nday1=count1
#	offset_day=offset*365
#	nday=count+2*offset_day
	
	origin <- paste(start_year-offset,"1","1",sep="/")
	if (yearly) {
		no_mean=TRUE
	} else {
		no_mean=FALSE
	}
	output <- splineInterpolateMonthlytoDaily(nday=nday,val=tval,origin=origin,first_row=offset_day+1,last_row=nday1,no_spline=no_spline,no_mean=no_mean)
	
	if (!yearly) {
		
		out <- output[(offset_day+1):nday1,]
		
	} else {
		
		origin_rescaling <- paste(start_year,"1","1",sep="/")
		out <- rescaling_monthly(data=output[(offset_day+1):nday1,],origin=origin_rescaling,val=val)
		
	}
	
	return(out)

}

