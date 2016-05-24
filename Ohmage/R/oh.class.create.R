#' Create a new class
#' @param class_urn class id
#' @param class_name class name
#' @param description class description
#' @param ... other arguments passed to ohmage
#' @export
oh.class.create <- function(class_urn, class_name, description="My new class.", ...){
	
	xhr <- oh.call("/class/create ", class_urn = class_urn, class_name = class_name, description=description, ...);		
	if(xhr$result == "success") {
		message("Class ", class_name, " created!");
	} else {
		return(xhr);
	}		
}
