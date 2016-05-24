#' phenology_MHmcmc_p generates set of parameters to be used with MHmcmc()
#' @title Generates set of parameters to be used with phenology_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a fit_phenology() fit
#' @param accept If TRUE, does not wait for use interaction
#' @description Interactive script used to generate set of parameters to be used with phenology_MHmcmc().\cr
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#'   	reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Gratiot"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters=3, las=1, xlim=c(-10, 300))
#' }
#' @export

phenology_MHmcmc_p<-function(result=stop("An output from fit_phenology() must be provided"), accept=FALSE) {

if (class(result)!="phenology") {
  stop("An output from fit_phenology() must be provided")
}

# d'abord je sors les paramètres à utiliser

par <- result$par

# "Peak"
Peak <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["Peak"]), 160, par["Peak"]))

# "Flat"
Flat <- c("dunif", 0, 50, 2, 0, 50, ifelse(is.na(par["Flat"]), 5, par["Flat"]))

# "Begin"
Begin <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["Begin"]), 160, par["Begin"]))

# "End"
End <- c("dunif", 0, 365, 2, 0, 365, ifelse(is.na(par["End"]), 220, par["End"]))

# "Length"
Length <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["Length"]), 100, par["Length"]))

# "LengthE"
LengthE <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["LengthE"]), 100, par["LengthE"]))

# "LengthB"
LengthB <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["LengthB"]), 100, par["LengthB"]))

# "Max"
Max <- c("dunif", 0, 200, 2, 0, 200, ifelse(is.na(par["Max"]), 100, par["Max"]))

# "PMin"
PMin <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMin"]), 2, par["PMin"]))

# "Min"
Min <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["Min"]), 1, par["Min"]))

# "PMinE"
PMinE <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMinE"]), 2, par["PMinE"]))

# "MinE"
MinE <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["MinE"]), 1, par["MinE"]))

# "PMinB"
PMinB <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["PMinB"]), 2, par["PMinB"]))

# "MinB"
MinB <- c("dunif", 0, 5, 2, 0, 5, ifelse(is.na(par["MinB"]), 1, par["MinB"]))

# "Phi"
Phi <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi"]), 10, par["Phi"]))

# "Delta"
Delta <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta"]), 10, par["Delta"]))

# "Alpha"
Alpha <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha"]), 10, par["Alpha"]))

# "Beta"
Beta <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta"]), 10, par["Beta"]))

# "Tau"
Tau <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau"]), 1, par["Tau"]))

# "Phi1"
Phi1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi1"]), 10, par["Phi1"]))

# "Delta1"
Delta1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta1"]), 10, par["Delta1"]))

# "Alpha1"
Alpha1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha1"]), 10, par["Alpha1"]))

# "Beta1"
Beta1 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta1"]), 10, par["Beta1"]))

# "Tau1"
Tau1 <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau1"]), 1, par["Tau1"]))

# "Phi2"
Phi2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Phi2"]), 10, par["Phi2"]))

# "Delta2"
Delta2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Delta2"]), 10, par["Delta2"]))

# "Alpha2"
Alpha2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Alpha2"]), 10, par["Alpha2"]))

# "Beta2"
Beta2 <- c("dunif", 0, 30, 2, 0, 30, ifelse(is.na(par["Beta2"]), 10, par["Beta2"]))

# "Tau2"
Tau2 <- c("dunif", 0, 2, 2, 0, 2, ifelse(is.na(par["Tau2"]), 1, par["Tau2"]))

# "Theta"
Theta <- c("dunif", 0, 10, 2, 0, 10, ifelse(is.na(par["Theta"]), 10, par["Theta"]))


priors <- list(Peak, Flat, Begin, End, Length, LengthE, LengthB, 
Length, PMin, Min, PMinE, MinE, PMinB, MinB, Phi, Delta, Alpha, 
Beta, Tau, Phi1, Delta1, Alpha1, Beta1, Tau1, Phi2, Delta2, 
Alpha2, Beta2, Tau2, Theta)

names(priors) <- c("Peak", "Flat", "Begin", "End", "Length", 
"LengthE", "LengthB", "Length", "PMin", "Min", "PMinE", "MinE", 
"PMinB", "MinB", "Phi", "Delta", "Alpha", "Beta", "Tau", "Phi1", 
"Delta1", "Alpha1", "Beta1", "Tau1", "Phi2", "Delta2", "Alpha2", 
"Beta2", "Tau2", "Theta")

for (i in 1:length(par)) {

		if (substr(names(par[i]), 1, 4)=="Max_") {
			priors <- c(priors, list(c("dunif", 0, 200, 2, 0, 200, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		
		if (substr(names(par[i]), 1, 4)=="Min_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		if (substr(names(par[i]), 1, 5)=="MinE_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		
		if (substr(names(par[i]), 1, 5)=="MinB_") {
			priors <- c(priors, list(c("dunif", 0, 5, 2, 0, 5, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
		if (substr(names(par[i]), 1, 5)=="Peak_") {
			priors <- c(priors, list(c("dunif", 0, 365, 2, 0, 365, par[i])))
			names(priors)[length(priors)] <- names(par[i])
		}
  if (substr(names(par[i]), 1, 5)=="Length_") {
    priors <- c(priors, list(c("dunif", 0, 100, 2, 0, 100, par[i])))
    names(priors)[length(priors)] <- names(par[i])
  }
}


prencours <- NULL

for (i in 1:length(par)) {
	prencours <- c(prencours, priors[[names(par)[i]]])
}



parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
rownames(parametersMCMC)<-names(par)
parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)

for (i in 2:7)
  parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])


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
	if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
	cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
	
	} else {
	
	if (density == "dnorm") {
	
	cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
	cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)
	
	} else {

	cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- as.numeric(f)
	cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- as.numeric(f)

	}
	}
	
	
	cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "SDProp"] <- as.numeric(f)
	cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Min"] <- as.numeric(f)
	cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Max"] <- as.numeric(f)
	cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Init"] <- as.numeric(f)
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
    warning(paste0("Initial value for parameter ", rownames(parameters)[i], " was out of range; It has been corrected. Check it.")) 
  }
}


}
