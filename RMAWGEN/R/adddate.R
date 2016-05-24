NULL
#'
#'  Inserts three columns (year,month,day) passing dates  to a matrix or to a dataframe
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#' @param data matrix of daily data 
#' @param origin character string containing the date of the first row of \code{data} as YYYY-MM-DD
#' 
#' @seealso \code{\link{findDate}}
#' @export 
#' @return a data frame with dates and \code{data} values
#' 


adddate <- function(data,origin="1961-1-1") {
	
	
	dates <- findDate(1:nrow(data),origin=origin,data.frame=TRUE,decimal=FALSE,character=FALSE)
	
	out <- cbind(dates,data)
	
	
	
	return(out)
}	
