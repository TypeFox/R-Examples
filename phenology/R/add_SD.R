#' add_SD adds SD for a fixed parameter.
#' @title Add SD for a fixed parameter.
#' @author Marc Girondot
#' @return The parameters set with the new SD value
#' @param parametersfixed Set of fixed parameters
#' @param parameters Set of current parameters
#' @param SD Standard deviation value to be added
#' @param help If TRUE, an help is displayed
#' @description This function is used to add standard deviation for a fixed parameter.
#' @examples
#' library(phenology)
#' # Generate a set of fixed parameter: Flat and Min
#' pfixed<-c(Flat=0, Min=0)
#'	# Add SD for the Flat parameter
#' pfixed<-add_SD(parametersfixed=pfixed, parameters="Flat", SD=5)

#' @export


add_SD <-
function(parametersfixed=NULL, parameters=NULL, SD=NULL, help=FALSE) {
if(help) {
	cat("This function is used to add standard deviation for a fixed parameter.\n")
	cat("The syntax is parfixed<-add_SD(parametersfixed=NULL, parameters=name, SD=value)\n")

} else {
if (!is.null(parametersfixed) && !is.null(SD) && !is.null(parameters)) {
	c<-SD
	names(c)<-paste("sd#",parameters,sep="")
	return(c(parametersfixed, c))
}
}
}
