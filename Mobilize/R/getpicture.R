#' Wrapper for oh.image.read that adds a contenttype
#' @param ... arguments passed on to oh.image.read
#' @return file with the photo image.
#' @export
getpicture <- function(...){
	myimage <- oh.image.read(...);
	
	attr(myimage, "CONTENTTYPE") <- "image/jpeg";
	return(myimage);
}
