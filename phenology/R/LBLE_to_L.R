#' LBLE_to_L transforms a set of parameters from LengthB LengthE format to Length.
#' @title Transform a set of parameters from LengthB LengthE format to Length
#' @author Marc Girondot
#' @return Return the set of modified parameters
#' @param parameters Set of current parameters
#' @description This function is used to transform a set of parameters 
#' that uses LengthB and LengthE to a set of parameters uses Length.
#' @examples
#' \dontrun{
#' # Read a file with data
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", 
#' 		header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Change the parameters to Begin End format
#' parg1<-LBLE_to_L(parameters=parg)
#' # And change back to LengthB LengthE.
#' parg2<-L_to_LBLE(parameters=parg1)
#' @export

LBLE_to_L <-
function(parameters=stop("Set of parameters must be given")) {

if (!is.na(parameters["LengthB"]) & !is.na(parameters["LengthE"])) {
	
	parameters["Length"] <- (parameters["LengthB"]+parameters["LengthE"])/2
	parameters <- parameters[-which(names(parameters)=="LengthB")]
	parameters <- parameters[-which(names(parameters)=="LengthE")]

	return(parameters)
	}

}
