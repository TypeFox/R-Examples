#' Update a class
#' @param class_urn class id 
#' @param ... arguments passed to ohmage
#' @export
oh.class.update <- function(class_urn, ...){
	
	xhr <- oh.call("/class/update ", class_urn = class_urn, ...);		
	if(xhr$result == "success") {
		message("Class updated!");
	} else {
		return(xhr);
	}		
}

