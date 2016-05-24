#' @export
as.precintcon.monthly <- function(object) {
	
	if (is.element("precintcon.daily", class(object)))
		return(precintcon.monthly.aggregation(object))
	
	else if (is.element("data.frame", class(object))) {
		
		if (ncol(object) == 3) {
			
			if (is.element(class(as.vector(object[,3])), c("integer", "numeric")) && 
				is.element(class(as.vector(object[,2])), c("integer", "numeric")) && 
				is.element(class(as.vector(object[,1])), c("integer", "numeric"))) {
				
				class(object) <- c("data.frame", "precintcon.monthly")
				colnames(object) <- c("year", "month", "precipitation")
				return(object)
			
			}
		}
		
	} else
		stop("Invalid input object. It should be of type data.frame or precintcon.daily and has only numeric values.")
}