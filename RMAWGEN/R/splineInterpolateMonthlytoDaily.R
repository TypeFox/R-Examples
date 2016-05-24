
NULL

#' 
#' Interpolates monthly data to daily data using \code{\link{spline}}  and preserving monthly mean values
#' 
#' @param nday number of days on which the daily data is requested, e.g. number of days in one year 
#' @param val matrix containing monthly mean data
#' @param origin  date corresponding to the first row of the returned matrix
#' @param first_row row corresponding the first day of time interval where montlhy mean conservation is applied
#' @param last_row corresponding the last day of time interval where montlhy mean conservation is applied
#' @param no_spline logical value. If \code{TRUE} no spline interpolation is calculated and the daily value corresponds to the monthly average value. Default is \code{FALSE}.
#' @param no_mean logical value. Default is \code{FALSE}. If \code{TRUE} the function output is not rescaled in order to maintain observed mean monthly values. 
#'  @export 
#'  
#' @return a matrix or data frame with interpolated daily data 
#' 
#' @seealso \code{\link{spline}},\code{\link{splineInterpolateMonthlytoDailyforSeveralYears}}
#' 
#' @author Emanuele Cordano, Emanuele Eccel

splineInterpolateMonthlytoDaily <-
function(nday=365,val=as.matrix(cbind(1*(0.5:11.5)*nday/12,2*(0.5:11.5)*nday/12)),origin="1961-1-1",first_row=1,last_row=nday,no_spline=FALSE,no_mean=FALSE) {
	


	nmonth=nrow(val)

	
	frac <- as.double(nday/nmonth)
	x <- ((1:nmonth)-0.5)*frac
	x_out <- 1:nday
	
	i_out <- as.integer((x_out-1)/frac)+1
	
	
	
	
	output <- array(NA,c(nday,ncol(val)))
	
	sval <- as.matrix(array(NA,c(nday,ncol(val))))
	

	
	
	for (i in 1:ncol(val)) {
		out <- spline(x=x,y=val[,i],xout=x_out)
		
		
		out1 <- as.vector(out$y) # insert monthds correctly!!!
		out2 <- out1
		
		months <- months((0.5:11.5)*365/12,abbreviate=TRUE)

		for (m in 1:length(months)) {
			
			i_months <- extractmonths(data=1:length(out2),when=months[m],origin=origin)
			
			
			
			if (no_spline) {
				out2[i_months[i_months %in% first_row:last_row]] <-  val[m,i]
			} else if (!no_mean){
				out2[i_months[i_months %in% first_row:last_row]] <-  out1[i_months[i_months %in% first_row:last_row]]-mean(out1[i_months[i_months %in% first_row:last_row]])+val[m,i] 
			} else {
				out2[i_months[i_months %in% first_row:last_row]] <-  out1[i_months[i_months %in% first_row:last_row]]
			}
	 	}	
#		for (m in 1:nmonth) {
			
#			out2[i_months] <-  out1[i_months]-mean(out1[i_months])+val[m,i] 
	
#		}
#		
		
		output[,i] <- out2
	}
	
	
	
	
	
	return(output)
	
}

