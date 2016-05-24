#' fitRMU_MHmcmc_p generates set of parameters to be used with fitRMU_MHmcmc()
#' @title Generates set of parameters to be used with fitRMU_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a fitRMU() fit
#' @param accept If TRUE, does not wait for use interaction
#' @family Fill gaps in RMU
#' @description Interactive script used to generate set of parameters to be used with phenology_MHmcmc().\cr
#' @examples 
#' \dontrun{
#' library("phenology")
#' RMU.name.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
#'                                          "Galibi.Suriname", 
#'                                          "Irakumpapy.French.Guiana"), 
#'                                  se=c("se_Yalimapo.French.Guiana", 
#'                                       "se_Galibi.Suriname", 
#'                                       "se_Irakumpapy.French.Guiana"))
#' data.AtlanticW <- data.frame(Year=c(1990:2000), 
#'       Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
#'                                6542, 5678, 1243, NA, 1566, 1566),
#'       se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
#'                                  230, 129, 167, NA, 145, 20),
#'       Galibi.Suriname=c(276, 275, 290, NA, 267, 
#'                        542, 678, NA, 243, 156, 123),
#'       se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
#'                            4.3, 2.3, NA, 10.3, 10.1, 8.9),
#'       Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
#'                                3542, 2678, 243, NA, 566, 566),
#'       se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
#'                                  130, 29, 67, NA, 15, 20))
#'                            
#' cst <- fitRMU(data=data.AtlanticW, RMU.name=RMU.name.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#' pMCMC <- fitRMU_MHmcmc_p(result=cst, accept=TRUE)
#' }
#' @export

fitRMU_MHmcmc_p <- function(result=stop("An output from fitRMU() must be provided"), 
	accept=FALSE) {

if (class(result)!="fitRMU") {
  stop("An output from fitRMU() must be provided")
}


# rownames(parametersMCMC)<-names(par)

parametersMCMC <- data.frame(Density=character(), Prior1=numeric(), Prior2=numeric(), 
SDProp=numeric(), Min=numeric(), Max=numeric(), Init=numeric(), stringsAsFactors = FALSE)

for (indice.par in seq_along(result$par)) {
  par <- result$par[indice.par]
  SE <- result$SE[indice.par]
  if (is.na(SE)) SE <- par/10
  nm <- names(par)
  if (nm=="r") parametersMCMC_ec <- data.frame(Density="dnorm", Prior1=unname(par), Prior2=unname(SE), 
                                            SDProp=2, Min=min(-1, unname(par)-1), 
                                            Max=max(1, unname(par)+1), Init=unname(par), 
                                            row.names = nm, stringsAsFactors = FALSE)
  if (substr(nm, 1, 2)=="T_") parametersMCMC_ec <- data.frame(Density="dnorm", 
                                              Prior1=unname(par), Prior2=unname(SE), 
                                               SDProp=2, Min=0.01, Max=max(10000, unname(par)+1000), Init=unname(par), 
                                      row.names = nm, stringsAsFactors = FALSE)
  
  if (substr(nm, 1, 3)=="SD_") parametersMCMC_ec <- data.frame(Density="dnorm", 
                                                              Prior1=unname(par), Prior2=unname(SE), 
                                                              SDProp=2, Min=0.01, Max=max(200, unname(par)+100), Init=unname(par), 
                                                              row.names = nm, stringsAsFactors = FALSE)
  if (substr(nm, 1, 1)=="a" & substr(nm, 3, 3)=="_") parametersMCMC_ec <- data.frame(Density="dnorm", 
                                                              Prior1=unname(par), Prior2=unname(SE), 
                                                              SDProp=2, Min=min(-10, unname(par)-2), 
                                                              Max=max(10, unname(par)+2), Init=unname(par), 
                                                              row.names = nm, stringsAsFactors = FALSE)
  parametersMCMC <- rbind(parametersMCMC, parametersMCMC_ec)
}

parametersMCMC[, "Prior2"] <- abs(parametersMCMC[, "Prior2"])

parameters <- parametersMCMC

if (accept) {
	return(parameters)
} else {

	repeat {

cat("Proposition:\n")
print(parameters)
cat("Name of the parameter to change or Enter to quit:\n")
f<-scan(nmax=1, quiet=TRUE, what=character())

if (length(f)==0) f <- "q"

if (f=="q") {
	return(parameters)
	
} else {

	variable <- which(f==names(par))
	if (length(variable)==0) {
	cat("The parameter does not exist:\n")
	} else {
	print(variable)
	cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))

	cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
	density<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
	
	if (density == "dunif") {
	
	cat(paste("Distribution of the prior, Minimum (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {
	
	if (density == "dnorm") {
	
	cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {

	cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f

	}
	}
	
	
	cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "SDProp"] <- f
	cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Min"] <- f
	cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Max"] <- f
	cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Init"] <- f
	}

}

}

}

for (i in 1:nrow(parameters)) {
  if (parameters[i, "Density"]=="dunif") {
    mn <- max(as.numeric(parameters[i, "Prior1"]), as.numeric(parameters[i, "Min"]))    
  } else {
    mn <- as.numeric(parameters[i, "Min"])
    mx <- as.numeric(parameters[i, "Max"])
  }  
  if (findInterval(as.numeric(parameters[i, "Init"]), c(mn, mx)) != 1) {
    parameters[i, "Init"] <- as.character(mn+(mx-mn)/2)
    warning(paste("Initial value for parameter ", rownames(parameters)[i], " was out of range; It is corrected. Check it.")) 
  }
}


}
