#' This method builds an object that carries necessary configuration values. 
#' The resulting object is a list, which you can modify outside of this function.
#' Currently contains tsHost, tsPort, openField,closeField, highField, lowField
#' and volField.
#' 
#' @return This function returns a plain list with configuration settings.
aqInit <- function(){
	# used to set default parameters. 
	# all these can be overriden. 

	ret = list()
	ret$tsHost = "127.0.0.1"
	ret$tsPort = 44444

  ret$openField = "OPEN"
	ret$closeField = "CLOSE" 
	ret$highField = "HIGH"
	ret$lowField = "LOW"
	ret$volField = "VOLUME"

	return(ret)

}


