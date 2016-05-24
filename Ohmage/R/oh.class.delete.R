#' Delete a class from the server
#' @param class_urn class id
#' @param ... other args passed to ohmage 
#' @export
oh.class.delete <- function(class_urn, ...){
	xhr <- oh.call("/class/delete",class_urn=class_urn, ...);	
	message("Class deleted!");
}
