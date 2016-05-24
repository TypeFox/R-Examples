residuals.sma <- function(object, ...){
	return(fitted(object, type = "residuals",...))
}