NULL

#' 
#' @export 
#' @method generate PrecipitationAmountModel
#' @S3method generate PrecipitationAmountModel
#' @aliases generate generate.PrecipitationAmountModel
#' @rdname generate
#'
generate.PrecipitationAmountModel <- function(x,...) {
	
	
	out <- predict(x,precipitation.value.random.generation=TRUE,...)
	
	
	return(out)
	
}
	