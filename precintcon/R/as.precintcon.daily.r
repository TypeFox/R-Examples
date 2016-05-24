#' @export
as.precintcon.daily <- function(object, na.value = NA) {

	if (!is.element("data.frame", class(object)) || ncol(object) != 33)
		stop("Invalid object. It must be a \"data.frame\" and contain 33 numeric columns.")
	
	for (i in 3:33) {
		for (j in 1:nrow(object)) {
			if (is.na(as.Date(paste(object[j,1], object[j,2], i - 2, sep = "-"), "%Y-%m-%d"))) {
				a <- NA
				object[j,i] <- a  
			}
		}
	}
	
	colnames(object) <- c("year", "month", paste(rep("d", 31), 1:31, sep=""))
	
	class(object) <- c("data.frame", "precintcon.daily")
	
	return(object)
		
}
