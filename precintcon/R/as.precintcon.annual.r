#' @export
as.precintcon.annual <- function(object) {
	
	if (!is.element("precintcon.daily", class(object)) && !is.element("precintcon.monthly", class(object)))
		stop("input data should be either of type precintcon.monthly or precintcon.daily")
	
	if (is.element("precintcon.daily", class(object)))
		object <- as.precintcon.monthly(object)
	
	result <- aggregate(object[,3], by = list(object[,1]), FUN = sum)
	
	colnames(result) <- c("year", "precipitation")
	
	class(result) <- c("data.frame", "precintcon.annual")
	
	return(result)
}