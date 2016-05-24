NULL


#' 
#' This function adjusts the monthly mean to a daily weather dataset (e. g. spline-interpolated temperature)
#' 
#' @param data data frame of wheather variables) 
#' @param val monthly means returned by \code{\link{getMonthlyMean}}
#' @param origin character string containing the gregorian date of the first day of \code{data} 
#' @author  Emanuele Cordano
#'
#' 
#'  @export    
#'   
#'
#' @return   A data frame with data of \code{data} rescaled with \code{val} for each month
#'      
#' @note It uses \code{\link{months}} and  \code{\link{julian}}
#' @seealso \code{\link{extractdays}}
#' 
#' 






rescaling_monthly <- function(data,val,origin="1961-1-1"){
# the folloeing commented lines are DEPRACETED
	#function(data=array(1:ndim_max,dim=c(ndim_max,1)),val,ndim_max=100000,origin="1961-1-1"){
#	str(data)
#	print(class(data))
	out <- adddate(as.data.frame(data),origin=origin)
	
	years <-  unique(out$year) #as.numeric(as.character(years(as.chron(as.POSIXct(origin,tz="GMT")+1:ndata-1))))
	months <- months((0.5:11.5)*365/12,abbreviate=TRUE)
	
	
	for (m in 1:length(months)) {
	
		i_months <- extractmonths(data=1:nrow(out),when=months[m],origin=origin)
		temp0 <- out[i_months,]

		for (year in years) {
			
			
			temp <- temp0[temp0$year==year,-(1:3)]
			
			for (i in 1:ncol(temp)) {
				
				mean <- val[[as.character(year)]][m,i]
				
				vtemp <- temp[,i]-mean(temp[,i],na.rm=TRUE)+mean
			
			
			
				temp[,i] <- vtemp
				
			}
	
			temp0[temp0$year==year,-(1:3)] <- temp
			
			
			
		}
		out[i_months,] <- temp0
		
	}
	
	return(out[,-(1:3)])
	
}
	
