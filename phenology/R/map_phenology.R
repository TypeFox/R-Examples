#' map_phenology generates a likelihood map.
#' @title Generate a likelihood map varying Phi and Delta.
#' @author Marc Girondot
#' @return Display a likelihood map
#' @param data dataset generated with add_format
#' @param parametersfixed Set of fixed parameters
#' @param parametersfit Set of parameters to be fitted
#' @param Phi Phi values to be analyzed
#' @param Delta Delta value to be analyzed
#' @param method_incertitude 'combinatory' estimates likelihood of all combinations for nest numbers;\cr
#'                           'convolution' [default] uses the exact likelihood of the sum of negative binomial distribution.
#' @param infinite Number of iterations for dmnbinom() used for method_incertitude='convolution'
#' @param zero_counts Example c(TRUE, TRUE, FALSE) indicates whether the zeros have 
#'                    been recorded for each of these timeseries. Defaut is TRUE for all.
#' @param progressbar If FALSE, do not show the progress bar
#' @description This function generates a map of likelihood varying Phi and Delta.\cr
#' 	Parameters are the same than for the fit_phenology() function except for trace that is disabled.\cr
#' 	If Alpha, Beta or Tau are not indicated, Alpha and Tau are set to 0 and 1 and Beta is fitted.\cr
#' 	Only one set of Alpha, Beta, Tau, Phi and Delta are used for all timeseries present in data.\cr
#' 	Note that it is possible to fit or fixed Alpha[n], Beta[n], Tau[n], Phi[n] and Delta[n] with [n]=1 or 2 
#' 	and then it is possible to use this function to establish the likelihood map for a 
#' 	second or third sinusoids added to the global pattern.\cr
#' 	If Delta is not specified, it is estimated from Phi and the same precision as Phi is used.\cr
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
#' map_Gratiot<-map_phenology(data=data_Gratiot, Phi=seq(from=0.1, to=20, length.out=100), 
#' 		parametersfit=parg2, parametersfixed=pfixed)
#' }
#' data(map_Gratiot)
#' # Plot the map
#' plot(map_Gratiot, col=heat.colors(128))
#' # Plot the min(-Ln L) for Phi varying at any delta value
#' plot_phi(map=map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi equal to the value for maximum likelihood
#' plot_delta(map=map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi the nearest to 15
#' plot_delta(map=map_Gratiot, Phi=15)
#' @export

map_phenology <-
function(data=NULL, parametersfit=NULL, parametersfixed=NA, 
         Phi=seq(from=0.2,to=20, length.out=100), Delta=NULL, infinite=50, 
         method_incertitude="convolution", zero_counts=TRUE, progressbar=TRUE) {

#.phenology.env<- NULL
#rm(.phenology.env)
  
  method_incertitude <- tolower(method_incertitude)
  if (method_incertitude=="convolution") method_incertitude <- 1
  if (method_incertitude=="combinatory") method_incertitude <- 2

if (is.null(parametersfixed)) {parametersfixed<-NA}
if (is.null(parametersfit)) {parametersfit<-NA}

# je varie Phi et Delta

#create 2 vectors in form of numeric sequence, for Delta and Phi
Phivalue=Phi
if (is.null(Delta)) {

Deltavalue=seq(from=0, to=max(Phivalue)/2, length.out=length(Phivalue)+1)

} else {
	Deltavalue=Delta
}

LPhi<-length(Phivalue)
LDelta<-length(Deltavalue)

# SET MATRIX
matrix(data=NA, LPhi, LDelta) -> input

	if (length(zero_counts)==1) {zero_counts<-rep(zero_counts, length(data))}
	if (length(zero_counts)!=length(data)) {
		stop("zero_counts parameter must be TRUE (the zeros are used for all timeseries) or FALSE (the zeros are not used for all timeseries) or possess the same number of logical values than the number of series analyzed.")
	}


# si ni Alpha ni Beta ne sont à ajuster, je mets Beta
if (is.na(parametersfit["Alpha"]) && is.na(parametersfit["Beta"])) {
	if (all(is.na(parametersfit))) {
		parametersfit <- c(Beta=0)
	} else {
		parametersfit <- c(parametersfit, Beta=0)
	}
}

# si Beta est à la fois fixe et à ajuster, je le retire des fixes
if (!is.na(parametersfit["Beta"]) && !is.na(parametersfixed["Beta"])) {parametersfixed<-parametersfixed[!names(parametersfixed)=="Beta"]}

# je vérifie que Alpha, Beta et Tau apparaissent bien au moins une fois, sinon je les mets en fixe
xpar<-c(parametersfit, parametersfixed)
if (is.na(xpar["Alpha"])) {parametersfixed<-c(parametersfixed, Alpha=0)}
if (is.na(xpar["Beta"])) {parametersfixed<-c(parametersfixed, Beta=0)}
if (is.na(xpar["Tau"])) {parametersfixed<-c(parametersfixed, Tau=1)}

# Si Phi ou Delta sont indiqués en paramètres à ajuster, je les retire
parametersfit<-parametersfit[!names(parametersfit)=="Phi"]
parametersfit<-parametersfit[!names(parametersfit)=="Delta"]

# mais je les mets en paramètres fixes
if (is.na(parametersfixed["Phi"])) {parametersfixed<-c(parametersfixed, Phi=0)}
if (is.na(parametersfixed["Delta"])) {parametersfixed<-c(parametersfixed, Delta=0)}

	
if (progressbar) pb<-txtProgressBar(min=0, max=LDelta, style=3)

#parpre1<-parametersfit
parpre<-parametersfit

#FILLING MATRIX
for(j in 1:LDelta) {

XDelta<-Deltavalue[j]
for(i in 1:LPhi) {
  XPhi<-Phivalue[i]
  if (XDelta>=XPhi/2) {
    input[i,j]=NA
  } else {
  	parametersfixed["Delta"]<-XDelta
  	parametersfixed["Phi"]<-XPhi

#	assign("fixed", parametersfixed, envir=as.environment(.phenology.env))
    
    par<-parpre

    repeat {
    	    	
		resul<-optim(par, getFromNamespace(".Lnegbin", ns="phenology"), 
		             pt=list(data=data, fixed=parametersfixed, incertitude=method_incertitude, 
		                     zerocounts=zero_counts, out=TRUE, infinite=infinite), method="BFGS",control=list(trace=0, REPORT=1, maxit=500),hessian=FALSE)
		if (resul$convergence==0) break
		par <- resul$par
		# print("Convergence is not achieved. Optimization continues !")
	}
	
#	if (resul$value>378.717) {
#		print(resul$value)
#		print(parpre)
#		print(parametersfixed)
#	}
	
#  parpre<-resul$par
#  if (i==1) {parpre1<-parpre}
  input[i,j]<-resul$value
  }
#  Sys.sleep(0)
if (progressbar) setTxtProgressBar(pb, j)
}
}


Dv=as.vector(input)
pos=which.min(Dv)-1
j0=floor(pos/nrow(input))+1
i0=pos%%nrow(input)+1
print(paste("The minimum -Ln likelihood is ", input[i0, j0], sep=""))
print(paste("For Phi=",Phivalue[i0],sep=""))
print(paste("And Delta=",Deltavalue[j0],sep=""))

outputmap <- list(input=input, Phi=Phivalue, Delta=Deltavalue, Parametersfitted=names(parametersfit), 
Parametersfixed=parametersfixed, Data=names(data))

class(outputmap) <- "phenologymap"

return(outputmap)

}
