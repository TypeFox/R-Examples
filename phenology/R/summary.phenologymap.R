#' summary.phenologymap print information on a phenologymap object
#' @title Print information on a phenologymap object.
#' @author Marc Girondot
#' @return Return None
#' @param object A map generated with map_phenology.
#' @param ... Not used
#' @examples
#' library("phenology")
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
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed<-c(parg1, Alpha=0, Tau=1)
#' pfixed<-pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2<-c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map 
#' # [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' \dontrun{
#' map_Gratiot<-map_phenology(data=data_Gratiot, 
#' 		Phi=seq(from=0.1, to=20, length.out=100), 
#' 		parametersfit=parg2, parametersfixed=pfixed)
#' }
#' data(map_Gratiot)
#' # Print the information on a map
#' summary(map_Gratiot)
#' @method summary phenologymap
#' @export



summary.phenologymap <- function(object, ...) {

	cat("Fixed parameters:\n")
	for (i in 1:length(object$Parametersfixed)) {
		cat(paste(names(object$Parametersfixed[i]), "=", object$Parametersfixed[i], "\n", sep=""))
	}
	cat("Phi values:\n")
	cat(object$Phi)
	cat("\n")

	cat("Delta values:\n")
	cat(object$Delta)
	cat("\n")

}
