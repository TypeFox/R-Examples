#' GenerateTest Generate a data.frame that can be used as test value for searchR()
#' @title Generate a data.frame that can be used as test value for searchR()
#' @author Marc Girondot
#' @return A data.frame with size or mass at hatching for each nest
#' @param series Name of series or object from searchR()
#' @param size Size or mass at hatching. Will be recycled if necessary
#' @param previous Previous formated test data
#' @description Generate a data.frame that can be used as test value for searchR()
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' testsize1 <- GenerateTest(resultNest_4p)
#' testsize2 <- GenerateTest(series=resultNest_4p,  
#' 	size=c(Mean=39.3, SD=1.92))
#' }
#' @export


GenerateTest <-
function(series=stop("A result object or names of series must be provided"), size=NULL, previous=NULL) {

	if (is.null(size) & (class(series)!="NestsResult")) {
		stop("size or a result from searchR() must be provided")
	}
	
	if (class(series)=="NestsResult") {
		if (is.null(size)) testec <- series$test
		series <- series$data
	}
	
	if (class(series)=="Nests") {
		series <- names(series)
		series <- series[1:(length(series)-2)]
	}
	
	
	if (!is.null(size)) {
		if (is.na(size["Mean"]) | is.na(size["SD"]) | length(size["Mean"])!=length(size["Mean"])) {
			return("size must be a vector with same number of Mean and SD values")
		} else {
			mean <- rep(size["Mean"], length(series))[1:length(series)]
			sd <- rep(size["SD"], length(series))[1:length(series)]
			testec <- data.frame(Mean=mean, SD=sd, row.names=series)
		}
	}
	
	if (!is.null(previous)) {
		test <- rbind(previous, testec)
	} else {
		test <- testec
	}
		
	return(test)	

}
