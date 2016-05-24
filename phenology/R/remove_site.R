#' remove_site removes beach information from a set of parameters.
#' @title Removes site information from a set of parameters.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to remove the information of the site
#' from a set of parameters. It can be used to other timeseries after.
#' @examples
#' library(phenology)
#' # Read a file with data
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' }
#' data(result_Gratiot)
#' # Extract parameters form result
#' parg<-extract_result(result_Gratiot)
#' # Remove site information
#' parg1<-remove_site(parg)
#' @export



remove_site <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to remove the information of the site\n")
	cat("from a set of parameters. It can be used to other timeseries after.\n")

} else {
if (!is.null(parameters)) {
	for(i in 1:length(parameters)) {
		names(parameters)[i]<-strsplit(names(parameters[i]), "_")[[1]][1]
	}
	return(parameters)
}
}
}
