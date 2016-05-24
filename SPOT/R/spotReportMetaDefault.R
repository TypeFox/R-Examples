###################################################################################################
#'  Default Report for Meta Runs
#' 
#' Function to generate a simple report for meta runs.
#'
#' This function draws a scatterplot matrix (based on car), printing it to screen or pdf. 
#' If the \code{report.io.pdf} setting is TRUE
#' the graphic is printed to a pdf file (usually named like your .conf file, and placed in the same folder)
#' if \code{report.io.screen} is set TRUE the graphic is printed to the screen. Both can be FALSE or TRUE
#' at the same time. If the user does not specify those values, the defaults will be used as shown in
#' \code{\link{spotGetOptions}}, which means there will be only screen output, and no pdf.
#' 
#' @param spotConfig the configuration list of all spot parameters
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepReport}} 
#' @export
###################################################################################################
spotReportMetaDefault <- function(spotConfig) {		
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotReportMetaDefault")
	fbs.df <- read.table(spotConfig$io.fbsFileName			
			 , header = TRUE
			 , as.is=TRUE
	);
	print(summary(fbs.df))
	spotConfig
}	

