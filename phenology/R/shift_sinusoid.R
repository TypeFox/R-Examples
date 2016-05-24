#' shift_sinusoid shift sinusoid information.
#' @title Shift sinusoid information.
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters set of parameters
#' @param from The number of series to change
#' @param to The number of series to change
#' @param help If TRUE, an help is displayed
#' @description This function is used to shift sinusoid parameters from '', '1' or '2'.
#' @examples
#' # Read a file with data
#' library("phenology")
#' \dontrun{
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' }
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Fix parameter FLat to 0
#' pfixed=c(Flat=0)
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=pfixed)
#' # Fit is done
#' \dontrun{
#' result_Gratiot_Flat<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=pfixed, trace=1)
#' }
#' data(result_Gratiot_Flat)
#' parg<-extract_result(result_Gratiot_Flat)
#' # Add data for one sinusoid superimposed 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' parg<-c(parg, Alpha=0.5, Beta=0.8, Delta=3, Phi=15)
#' # Tau is fixed to 1
#' pfixed=c(Flat=0, Tau=1)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot1<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=pfixed, trace=1)
#' # Plot the phenology
#' output1<-plot(result_Gratiot1, moon=TRUE)
#' #' }
#' data(result_Gratiot1)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot1)
#' # Shift sunusoid information to the '1'
#' parg2<-shift_sinusoid(parameters=parg1, from="", to="1")
#' # Tau is fixed to 1
#' pfixed=c(Flat=0, Tau1=1, Tau=1)
#' # Add data for another sinusoid superimposed 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' parg<-c(parg2, Alpha=0.5, Beta=0.8, Delta=3, Phi=10)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot2<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=pfixed, trace=1)
#' # Plot the phenology
#' output2<-plot(result_Gratiot2, moon=TRUE)
#' }
#' data(result_Gratiot2)
#' @export


shift_sinusoid <-
function(parameters=NULL, from="", to="1", help=FALSE) {
if(help) {
	cat("This function is used to roll sinusoid parameters from '', '1' or '2'.\n")
	cat("The syntax is par<-shift_sinusoid(parameters=value, from='1', to='2')\n")
	cat("or\n")
	cat("par<-shift_sinusoid(parameters=value, from='', to='1')\n")

} else {
if (!is.null(parameters) && (from!=to)) {

	level<-from
	varfrom<-c(paste("Phi", level, sep=""), paste("Delta", level, sep=""),paste("Alpha", level, sep=""),paste("Beta", level, sep=""),paste("Tau", level, sep=""))
	level<-to
	varto<-c(paste("Phi", level, sep=""), paste("Delta", level, sep=""),paste("Alpha", level, sep=""),paste("Beta", level, sep=""),paste("Tau", level, sep=""))
	
	for (i in 1:5) {
		if (length(names(parameters)[names(parameters)==varfrom[i]])==1) {
			names(parameters)[names(parameters)==varfrom[i]]<-varto[i]
		}
	}
	return(parameters)

}
}
}
