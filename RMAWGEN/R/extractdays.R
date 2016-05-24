NULL


#' 
#' Extracts the rows of a matrix corresponding to the requested days (expressed as dates YYYY-MM-DD) given the  date (origin) of the first row 
#' 
#' @param data an input data matrix where each row corresponds to a daily record
#' @param when desired dates for which the data are requested
#' @param origin date corresponding to the first row of \code{data}
#' @param nday (optional) number of days since \code{when} to extract the data
#' @param ndim_max maximum (integer) number of rows in \code{data} where to find \code{when}. 
#' Default is 100000 and works if \code{data} is missing.
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'   
#' @export 
#' @note It uses \code{\link{julian}}
#'  
#'
#'        
#' 
#' @return  a matrix containing the requested rows



extractdays <-
function(data=array(1:ndim_max,dim=c(ndim_max,1)),ndim_max=100000,when="1990-1-1",origin="1961-1-1",nday=1){

	k <- as.integer(julian(as.POSIXct(when,tz="GMT"),origin=as.POSIXct(origin,tz="GMT")))
	
	
	if (min(k,na.rm=TRUE)<0) print("Error in extractdays: when parameter date is after the time origin")
	
	
	if (nday>1) {
		return(data[(k+1):(k+nday),])
	} 
	return(data[k+1,])
	
}

