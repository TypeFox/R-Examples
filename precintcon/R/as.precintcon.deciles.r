#' @export
as.precintcon.deciles <- function(object) {
	
	if (is.element("precintcon.deciles", class(object))) {
		return(object)
	} else {
		if (is.element("precintcon.daily", class(object)) || is.element("precintcon.monthly", class(object))) {
			object <- precintcon.deciles.analysis(object)
		} else if (!is.element("precintcon.deciles", class(object))){
			stop("invalid input data. It should be either of type precintcon.daily, precintcon.month, or precintcon.deciles.")
		}
		return(object)
	}
}
