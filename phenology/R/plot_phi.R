#' plot_phi plots the best likelihood for fixed Phi value.
#' @title Plot the best likelihood for fixed Phi value.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology
#' @param help If TRUE, an help is displayed
#' @description The function "plot_phi" plots the best likelihood for each Phi value.
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
#' # Generate a likelihood map [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' \dontrun{
#' map_Gratiot<-map_phenology(data=data_Gratiot, 
#' 		Phi=seq(from=0.1, to=20, length.out=100), 
#' 		parametersfit=parg2, parametersfixed=pfixed)
#' }
#' data(map_Gratiot)
#' # Plot the min(-Ln L) for Phi varying at any delta value
#' plot_phi(map=map_Gratiot)
#' @export


plot_phi <-
function(map=NULL, help=FALSE) {

if(help || is.null(map)) {
	cat("This function plots a likelihood lineplot obtained after map_phenology.\n")
	cat("The syntaxe is:\n")
	cat("plot_phi(map=mapx, pdf=TRUE, pdfname='NameMap.pdf')\n")

} else {

effetphi<-NULL
for(i in 1:length(map$Phi)) effetphi<-c(effetphi, min(map$input[i,], na.rm=TRUE))
plot(map$Phi, effetphi, type="l", xlab="Phi", ylab="-Ln L", bty="n", main=map$Data)


}
}
