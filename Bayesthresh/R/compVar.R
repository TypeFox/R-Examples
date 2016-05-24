# Extract the variance components of model

compVar <- function(object)
{
	if(!inherits(object, "Bayesthresh"))
		stop("Use an object of class Bayesthresh")
	return(object$compVar)
}
