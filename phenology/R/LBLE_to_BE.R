#' LBLE_to_BE transforms a set of parameters from LengthB LengthE to Begin End.
#' @title Transform a set of parameters from LengthB LengthE to Begin End.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of current parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to transform a set of parameters 
#' that uses LengthB, Peak and LengthE to a set of parameters 
#' that uses Begin, Peak and End.
#' @examples 
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Change the parameters to Begin End format
#' parg1<-LBLE_to_BE(parameters=parg)
#' # And change back to LengthB LengthE
#' parg2<-BE_to_LBLE(parameters=parg1)
#' @export


LBLE_to_BE <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to transform a set of parameters\n")
	cat("that uses LengthB, Peak and LengthE to a set of parameters\n")
	cat("that uses Begin, Peak and End.\n")

} else {

if (!is.null(parameters)) {
	pk<-parameters["Peak"]
	lb<-parameters["LengthB"]
	le<-parameters["LengthE"]

	bg<-pk-lb
	en<-pk+le
	
	parameters[names(parameters)=="LengthB"]<-bg
	names(parameters)[names(parameters)=="LengthB"]<-"Begin"
	
	parameters[names(parameters)=="LengthE"]<-en
	names(parameters)[names(parameters)=="LengthE"]<-"End"
	
	return(parameters)
	
}
}
}
