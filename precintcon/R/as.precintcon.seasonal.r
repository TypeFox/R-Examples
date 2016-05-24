#' @export
as.precintcon.seasonal <- function(object, hemisthere = c("n", "s")) {
	
  if (is.element("precintcon.seasonal", class(object)))
    return(object)
  
  object <- as.precintcon.monthly(object)
	
	result <- data.frame()
	
	season <- c("Spring", "Summer", "Autumn", "Winter")
	
	start <- which(object$month == 3)[1]
	
	for(i in seq(start, nrow(object), by = 3)) {
	  if (nrow(object) - i < 2) break;
	  
		result <- rbind(result, 
			data.frame(
				year          = object[i,1], 
				season        = season[abs((i - start) %/% 3 + ifelse(hemisthere == "n", 0, 2)) %% 4 + 1], 
				precipitation = sum(object[i:(i+2),3])))
	}

	class(result) <- c("data.frame", "precintcon.seasonal")

	return(result)
}