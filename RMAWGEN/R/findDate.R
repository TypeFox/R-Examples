
NULL


#' 
#' Finds the date  corresponding a row index of a matrix given the  date (origin) of the first row 
#' 
#' @param k integer or decimal value corresponding to number of days since \code{origin}
#' @param origin origin date. See also \code{\link{extractdays}}
#' @param data.frame logical variable. If \code{TRUE} (default) the date is returned as data frame (like \code{data} in \code{\link{extractyears}}), otherwise it is returned as  character or POSIXct. 
#' @param decimal logical variable. If \code{FALSE} (default) k is integer and starts from 1, otherwise is consider as the decimal julian day since \code{origin} (deprecated)
#' @param character logical variable. It is used if \code{data.frame} is \code{FALSE}, if it is \code{FALSE}, the date is returned as POSIXct, otherwise it is a character in the following form: \code{YYYY-MM-DD} 
#' 
#' @export
#' 
#'  @author  Emanuele Cordano, Emanuele Eccel
#'    
#' @import date 
#'
#'         
#' 
#' @note It uses functions of  \code{\link{time}} package. It works like an  inverse functions of \code{\link{extractdays}}. 
#' If \code{k} is a vector, the function returns several dates for each element of \code{k}
#' 
#' 
#' 
#' @seealso \code{\link{date.mdy}},\code{\link{extractdays}}
#' 
#' 
#' 
#' @return   the date(s) corresponding to \code{k} under different formats
#' 
#' @examples 
#' findDate <- findDate(100,origin="1961-1-1",data.frame=FALSE,character=TRUE)
#' 
#' 





findDate <-
function(k,origin="1961-1-1",data.frame=TRUE,decimal=FALSE,character=FALSE) {

	
	origin0 <- as.POSIXct(paste(date.mdy(0)$year,date.mdy(0)$month,date.mdy(0)$day,sep="-"),tz="GMT")
	
	
	
	k_origin=as.integer(julian(as.POSIXct(origin,tz="GMT"),origin=origin0))
	
	if (decimal) {
		out <- date.mdy(k+k_origin)
	} else {
		out <- date.mdy(k+k_origin-1)
	}
	if (data.frame) return(out)
	
	
	outdata <- paste(out$year,out$month,out$day,sep="-")
	
	if (!character) outdata <- as.POSIXct(paste(out$year,out$month,out$day,sep="-"),tz="GMT")
	
	if (character) outdata[is.na(out$year)]=array(NA,length(outdata[is.na(out$year)]))
	
	return(outdata)
}

