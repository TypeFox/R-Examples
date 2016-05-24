#' extract_result get the fitted parameters from a result object.
#' @title Extract the set of parameters from a result object.
#' @author Marc Girondot
#' @return Return the set of fitted parameters
#' @param result A result file
#' @param help If TRUE, an help is displayed
#' @description The function "extract_result" permits to extract the set of parameters from a result object obtained after fit_phenology.
#' @examples
#' library(phenology)
#' \dontrun{
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", \cr
#' 		header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, 
#' 		parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' }
#' @export


extract_result <-
function(result=NULL, help=FALSE) {
if(help) {
	cat("This function is used to get the set of parameters\n")
	cat("from a result object obtained after fit_phenology.\n")

} else {
	return(result$par)
}
}
