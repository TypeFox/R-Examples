#' Read info about a class
#' @param class_urn_list list of class id's
#' @param ... other arguments passed to ohmage
#' @export
oh.class.read <- function(class_urn_list="urn:sys:andwellness", ...){
	xhr <- oh.call("/class/read", class_urn_list=class_urn_list, ...);	
	return(xhr);	
}
