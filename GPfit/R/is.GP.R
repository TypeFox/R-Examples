## Checks whether or not the object in question is a GP object
##
## May 8th, 2012

is.GP <- function(object){
if (inherits(object, "GP"))
{
	return(TRUE)
} else {
	return(FALSE)
}
}