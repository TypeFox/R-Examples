################################################################################
# "Working with dynamic models for agriculture"
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The Magarey model
#' @description Generic model of infection for foliar diseases caused by fungi (from Magarey et al.,2005).
#' @param T : input variable. Either a scalar or a vector (for a weather series).
#' @param Tmin : parameter of minimal temperature for infection (degC)
#' @param Topt : parameter of optimal temperature for infection (degC)
#' @param Tmax : parameter of maximal temperature for infection (degC)
#' @param Wmin : parameter of minimal wetness duration for infection (hour)
#' @param Wmax : parameter of maximal wetness duration for infection (hour)
#' @return Wetness duration (W, hour). Either a scalar or a vector depending on T. 
#' @export
#' @examples plot(1:35, magarey.model (1:35,7, 18, 30, 10, 42), type="l", xlab="T", ylab="W")
magarey.model<-function(T, Tmin, Topt, Tmax, Wmin, Wmax){
    fT<-((Tmax-T)/(Tmax-Topt))*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))
	W <- Wmin/fT
	W[W>Wmax | T<Tmin | T>Tmax]<-Wmax 
	return(W)
}
#' @title The Magarey model, taking a vector of parameters as argument
#' @description Generic model of infection for foliar diseases caused by fungi (from Magarey et al.,2005).
#' @param param : parameters 
#' @param T : input variable. Either a scalar or a vector (for a weather series). 
#' @return W : Wetness duration (hour). Either a scalar or a vector depending on T.  
#' @export
magarey.model2<-function(T, param){
    return(magarey.model(T, param["Tmin"], param["Topt"], param["Tmax"], param["Wmin"], param["Wmax"]))
}


################################################################################
#' @title Define values of the parameters for the Magarey model
#' @param species : name of a species. By default, value for an "unkown" species are given. Other possibility are "G.citricarpa" or "Kaki.fungus"
#' @return matrix with parameter values (nominal, binf, bsup)
#' @examples magarey.define.param(species="G.citricarpa")
#' magarey.define.param(species="Kaki.fungus")
#' @export
magarey.define.param <- function(species="unkown")
{
# nominal, binf, bsup
# Tmin  : the minimal temperature for infection (°C)
# Topt  : the optimal temperature for infection (°C)
# Tmax  : the maximal temperature for infection (°C)
# Wmin : minimal wetness duration for infection (hour)
# Wmax : maximal wetness duration for infection (hour)

if (species=="unkown"){
Tmin = c(NA, 0.8, 1.2)
Topt = c(NA, 20, 30)
Tmax = c(NA, 24, 36)
Wmin = c(NA, 38.4, 57.6)
Wmax = c(NA, 115.2, 172.8)
print("See help for values for other species.")
}
# pycnidiospores of the fungus Guignardia citricarpa Kiely by EFSA (2008).
if (species=="G.citricarpa"){
Tmin = c(NA, 10, 15)
Topt = c(NA, 25, 30)
Tmax = c(NA, 32, 35)
Wmin = c(NA, 12, 14)
Wmax = c(NA, 35, 48)}
# Kaki fungus.
if (species=="Kaki.fungus"){
Tmin = c(7, 2, 13)
Topt = c(18, 14, 26)
Tmax = c(30, 27, 35)
Wmin = c(10, 5, 17)
Wmax = c(42, 18, 90)}
print(paste("parameter values for species : ",species))
param<-data.frame(Tmin,Topt,Tmax, Wmin, Wmax)
row.names(param)<-c("nominal","binf","bsup")
return(as.matrix(param))
}
################################################################################
#' @title Wrapper function to run the Magarey model multiple times (for multiple sets of inputs)
#' @param X : parameter matrix
#' @param T : input variable, temperature
#' @param all : if you want a matrix combining X and output
#' @return a table with wetness duration (W) for each parameter vector
#' @export
#' @description Example magarey.simule(magarey.define.param(),15)
magarey.simule <- function(X, T,  all=FALSE){
Y <- apply(X,1,function(param) magarey.model2(T, param))
if(all) Y = cbind(X,W = Y)
return(as.matrix(Y))
}
################################################################################
# End of file
