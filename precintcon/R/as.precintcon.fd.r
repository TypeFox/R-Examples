#' @export
as.precintcon.fd <- function(object, interval=1) {
	
	if (is.element("precintcon.fd", class(object)))
		return(object)
	
	else if (is.element("precintcon.classified", class(object)))
		return(precintcon.fd(object))
	
	else if (is.element("precintcon.daily", class(object)))
		return(precintcon.fd(precintcon.classification(object, interval)))
	
	else
		stop("invalid input data in as.")
}