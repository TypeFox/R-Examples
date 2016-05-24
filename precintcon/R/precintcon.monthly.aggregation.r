#' @noRd
#' @name precintcon.monthly.aggregation
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.monthly.aggregation 
#' @title Daily to monthly precipitation serie 
#' @description Transformes a daily precipitation serie into a monthly one. 
#' @usage precintcon.monthly.aggregation(object) 
#' @param object a daily precipitation serie.
#' @return A monthly precipitation serie calculated from a daily one. 
#' @seealso \code{\link{as.precintcon.monthly}} 
#' @keywords precipitation monthly granularity 
precintcon.monthly.aggregation <- function(object) {
	
	if ((is.element("precintcon.daily", class(object)))) {
		
		sum <- rowSums(object[3:33], na.rm=TRUE)
		
		data <- data.frame(object[1], object[2], sum)
		
		colnames(data) <- c("year", "month", "precipitation")
		
		class(data) <- c("data.frame", "precintcon.monthly")
		
		return(data)
		
	} else if ((is.element("precintcon.monthly", class(object)))) {
		
		return(object)
		
	} else
		stop("Invalid data. Please, check your input object.")
}