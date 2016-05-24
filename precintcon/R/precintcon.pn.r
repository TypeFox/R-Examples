#' @noRd
#' @name precintcon.pn
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.pn 
#' @title precintcon.pn 
#' @description precintcon.pn 
#' @usage precintcon.pn(x, interval = 30, scale = "a") 
#' @param x is a daily or monthly precipitation serie.
#' @param interval .
#' @param scale .
#' @return The PN index. 
#' @seealso \code{\link{precintcon.pn.analysis}} 
#' @keywords precipitation percent of normal 
precintcon.pn <- function(x, interval=30, scale="a") {
	
	if (!is.element("precintcon.daily", class(x)) && !is.element("precintcon.monthly", class(x)))
		stop("input data should be either of type precintcon.daily or precintcon.monthly")
	
	if (scale == "w") {
		
		stop("this method is not implemented for weekly scale")
		
		if (scale == "w" && !is.element("precintcon.daily", class(x)))
			stop("your dataset should be daily when using weekly scale")
		
	} else if (scale == "m")
		x <- as.precintcon.monthly(x)
	
	else if (scale == "s")
		x <- as.precintcon.seasonal(x)
	
	else if (scale == "a")
		x <- as.precintcon.annual(x)
	
	else
		x <- as.precintcon.decade(x)
	
	result <- data.frame()
	
	for (i in (interval):nrow(x)) {
		pmean <- mean(x[(i-interval+1):i,ncol(x)])
		pcurr <- x[i,ncol(x)]
		
		result <- rbind(result, cbind(if (scale != "a" && scale != "d") x[i,1:2] else x[i,1], pn=pcurr/pmean*100))
	}
	
	if (scale == "w") {
		colnames(result) <- c("year", "week", "pn")
	} else if (scale == "m") {
		colnames(result) <- c("year", "month", "pn")
	} else if (scale == "s") {
		colnames(result) <- c("year", "season", "pn")
	} else if (scale == "a") {
		colnames(result) <- c("year", "pn")
	} else if (scale == "d") {
		colnames(result) <- c("decade", "pn")
	}
	
	class(result) <- c("data.frame", "precintcon.pn")
	
	return(result)
}