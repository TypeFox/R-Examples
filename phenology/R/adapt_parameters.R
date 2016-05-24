#' adapt_parameters get the fitted parameters from a result object.
#' @title Extract the parameters from a set of parameters to be used with another dataset.
#' @author Marc Girondot
#' @return Return the set of parameters
#' @param parameters A set of parameters
#' @param data A dataset of counts
#' @description The function "adapt_parameters" extracts the set of parameters to be used with a subset of data. All the uncessary parameters are removed. It can be used when a set of beaches are fitted first and after only one of these beaches is fitted again.
#' @examples
#' library(phenology)
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=refdate, format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Add unnecessary parameters to parg
#' parg <- c(parg, Max_dummybeach=2, Peak_dummybeach=123)
#' # Extract the fitted parameters
#' parg1<-adapt_parameters(data=data_Gratiot, parameters=parg)
#' @export


adapt_parameters <-
function(data=stop("Datasets is mandatory for this function"), parameters=stop("Set of parameters is mandatory for this function")) {

x <- names(parameters)
x <- gsub("^Max_", "", x)
x <- gsub("^Peak_", "", x)
x <- gsub("^Min_", "", x)
x <- gsub("^MinB_", "", x)
x <- gsub("^MinE_", "", x)

y <- names(data)
y <- c(y, "Peak", "Length", "LengthB", "LengthE", "Begin", "End", "Phi", "Delta", "Alpha", "Beta", "Tau", 
"PMin", "PMinB", "PMinE",
"Phi1", "Delta1", "Alpha1", "Beta1", "Tau1",
"Phi2", "Delta2", "Alpha2", "Beta2", "Tau2", "Theta")

parred <- NULL

for (i in 1:length(y)) {

#	if (any(x)==y[i]) {
		parred <- c(parred, parameters[x==y[i]])
#	}
}

return(parred)

}
