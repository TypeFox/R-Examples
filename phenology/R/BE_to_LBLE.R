#' BE_to_LBLE transforms a set of parameters from Begin End format to LengthB LengthE.
#' @title Transform a set of parameters from Begin End to LengthB LengthE.
#' @author Marc Girondot
#' @return Return the set of modified parameters
#' @param parameters Set of current parameters
#' @param help If TRUE, an help is displayed
#' @description This function is used to transform a set of parameters 
#' that uses Begin, Peak and End to a set of parameters 
#' that uses LengthB, Peak and LengthE.
#' @examples
#' \dontrun{
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", 
#' 		header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Change the parameters to Begin End format
#' parg1<-LBLE_to_BE(parameters=parg)
#' # And change back to LengthB LengthE.
#' parg2<-BE_to_LBLE(parameters=parg1)
#' @export

BE_to_LBLE <-
function(parameters=NULL, help=FALSE) {
if(help) {
	cat("This function is used to transform a set of parameters\n")
	cat("that uses Begin, Peak and End to a set of parameters\n")
	cat("that uses LengthB, Peak and LengthE.\n")

} else {

if (!is.null(parameters)) {
	pk<-parameters["Peak"]
	bg<-parameters["Begin"]
	en<-parameters["End"]

	lb<-pk-bg
	le<-en-pk
	
	parameters[names(parameters)=="Begin"]<-lb
	names(parameters)[names(parameters)=="Begin"]<-"LengthB"
	
	parameters[names(parameters)=="End"]<-le
	names(parameters)[names(parameters)=="End"]<-"LengthE"
	
	return(parameters)
	
}
}
}
