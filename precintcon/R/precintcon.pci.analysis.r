#' @export 
precintcon.pci.analysis <- function(object) {
	
	if (is.element("precintcon.daily", class(object)))
		object <- as.precintcon.monthly(object);
	
	if (is.element("precintcon.monthly", class(object))) {
		
		p1 <- aggregate((object[,3])^2, by = list(object[,1]), FUN=sum)[2]
		
 		p2 <- aggregate(object[,3], by = list(object[,1]), FUN=sum)[2]
		
		data <- data.frame(year=unique(object[,1]), pci=(100 * (p1/(p2^2))))

		colnames(data) <- c("year", "pci")
		
		class(data) <- c("data.frame", "precintcon.pci")
		
		return(data)
		
	} else
		stop("Invalid data. Please, check your input object.")
}
