###################################################################################################
#' Earth Report
#' 
#' Generates a report for results of a SPOT run with spotPredictEarth
#'
#' @param spotConfig the configuration list of all spot parameters
#' @return list spotConfig with changed values
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepReport}} \code{\link{spotPredictEarth}} 
#' @export
###################################################################################################
spotReportEarth <- function(spotConfig) {	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotReportEarth")
	spotWriteLines(spotConfig$io.verbosity,2,"  ........................")
	spotWriteLines(spotConfig$io.verbosity,2,"Statistics:")
	spotPrint(spotConfig$io.verbosity,1,spotConfig$seq.modelFit)
	spotWriteLines(spotConfig$io.verbosity,2,"  ........................")
	spotWriteLines(spotConfig$io.verbosity,2,"Coefficients:")
	spotPrint(spotConfig$io.verbosity,1,spotConfig$seq.modelFit$coefficients)
	spotWriteLines(spotConfig$io.verbosity,2,"  ........................")
	if(class(spotConfig$seq.modelFit)=="earth"){
		spotWriteLines(spotConfig$io.verbosity,2,"Variable Importance Estimation Matrix:")
		spotPrint(spotConfig$io.verbosity,1,earth::evimp(spotConfig$seq.modelFit))
		spotWriteLines(spotConfig$io.verbosity,2,"  ........................")	
	}
	dev.new()
	plot(spotConfig$seq.modelFit)
	dev.new()
	plotmo::plotmo(spotConfig$seq.modelFit,all1=TRUE,all2=TRUE)
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotReportEarth")
	spotConfig		
}
