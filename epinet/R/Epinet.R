
# File Epinet.R

# FUNCTION MCMCcontrol
# Auxiliary function used to set MCMC control parameters to epinet()

MCMCcontrol <- function(nsamp, thinning, extrathinning = FALSE, burnin = 0, seed = floor(runif(1, 0, 2^30)), etapropsd)
{
	if (burnin > nsamp) stop("Error: Burn-in cannot be greater than number of samples requested.")
	if (burnin < 0) stop("Error: Burn-in cannot be negative.")
	out <- list(nsamp = nsamp, thinning = thinning, extrathinning = extrathinning, burnin = burnin, seed = seed, etapropsd = etapropsd)
	return(out)
}

# FUNCTION priorcontrol
# Auxiliary function used to set priors and hyperparameters for epinet()

priorcontrol <- function(bprior, tiprior, teprior, etaprior, kiprior, keprior, priordists = "gamma", 
	betapriordist = priordists, thetaipriordist = priordists, thetaepriordist = priordists, 
	etapriordist = "normal", kipriordist = priordists, kepriordist = priordists, 
	parentprobmult = 1)
{
	out <- list(bprior = bprior, tiprior = tiprior, teprior = teprior, etaprior = etaprior, 
				kiprior = kiprior, keprior = keprior, priordists = priordists, 
				betapriordist = betapriordist, thetaipriordist = thetaipriordist, 
				thetaepriordist = thetaepriordist, etapriordist = etapriordist, 
				kipriordist = kipriordist, kepriordist = kepriordist, parentprobmult = parentprobmult)
	return(out)
}

# FUNCTION epinet
# Formula interface for epibayesmcmc()
# Inputs: model formula, epidemic data, dyadic covariate matrix, MCMC controls, priors, and verbose flag
# Output: an object of class "epinet" containing posterior parameter samples

epinet <- function(formula, epidata, dyadiccovmat, mcmcinput = MCMCcontrol(), priors = priorcontrol(), verbose = TRUE)
{	
	saveseed <- mcmcinput$seed
	set.seed(saveseed)

	mf <- model.frame(formula = formula, data = as.data.frame(dyadiccovmat))
	mm <- model.matrix(attr(mf, "terms"), data = mf)
	newdyadiccovmat <-  cbind(dyadiccovmat[ , 1:2], mm)	
	
	ninf <- min(which(is.na(epidata[ , 5])),dim(epidata)[1] + 1) - 1
	mcmcinput$inferEtimes <- sum(is.na(epidata[1:ninf,3])) > 0	
	mcmcinput$inferItimes <- sum(is.na(epidata[1:ninf,4])) > 0	
	etapars <- dim(dyadiccovmat)[2] - 2
	if (length(priors$etapriordist) == 1) priors$etapriordist <- rep(priors$etapriordist, times = etapars)

	results <- epibayesmcmc(epidata = epidata, dyadiccovmat = newdyadiccovmat, nsamp = mcmcinput$nsamp,
				thinning = mcmcinput$thinning, bprior = priors$bprior, tiprior = priors$tiprior, 
				teprior = priors$teprior, etaprior = priors$etaprior, kiprior = priors$kiprior, 
				keprior = priors$keprior, etapropsd = mcmcinput$etapropsd, priordists = priors$priordists, 
				betapriordist = priors$betapriordist, thetaipriordist = priors$thetaipriordist, 
				thetaepriordist = priors$thetaepriordist, etapriordist = priors$etapriordist, 
				kipriordist = priors$kipriordist, kepriordist = priors$kepriordist, 
				extrathinning = mcmcinput$extrathinning, inferEtimes = mcmcinput$inferEtimes, 
				inferItimes = mcmcinput$inferItimes, parentprobmult = priors$parentprobmult, verbose = verbose,
				burnin = mcmcinput$burnin)	

	results$call <- match.call()
	results$formula <- formula
	results$mcmcinfo <- mcmcinput
	class(results) <- "epinet"	
	results		
}

# FUNCTION print.epinet
# Print method for class epinet

print.epinet <- function(x, ...)
{
	cat("Call: \n")
	print(x$call)
	cat("Network parameters: \n")
	for (i in 1:dim(x$eta)[2])
		cat("Eta", i, ": ", colnames(x$eta)[i], "\n")
	cat("\nMCMC chain: \n")	
	cat("Network and epidemic parameters:", x$mcmcinfo$nsamp, "iterations, thinned every", x$mcmcinfo$thinning, "iterations: ", length(x$llkd)[1], "posterior samples returned \n")
	if (x$mcmcinfo$inferEtimes) 
		{ if(!is.null(x$exptimes)) cat("Exposure times:", x$mcmcinfo$nsamp, "iterations, thinned every", x$mcmcinfo$thinning * x$mcmcinfo$extrathinning, "iterations: ",dim(x$exptimes)[2], "posterior samples returned \n") else cat("Exposure times inferred, but posterior samples not returned\n") }
	else cat("Exposure times provided, not inferred \n") 
	if (x$mcmcinfo$inferItimes) 
		{ if(!is.null(x$inftimes)) cat("Infection times:", x$mcmcinfo$nsamp, "iterations, thinned every", x$mcmcinfo$thinning * x$mcmcinfo$extrathinning, "iterations: ",dim(x$inftimes)[2], "posterior samples returned \n") else cat("Infection times inferred, but posterior samples not returned\n") }
	else cat("Infection times provided, not inferred \n") 
	if(!is.null(x$transtree)) cat("Transmission tree:", x$mcmcinfo$nsamp, "iterations, thinned every", x$mcmcinfo$thinning * x$mcmcinfo$extrathinning, "iterations: ",dim(x$transtree)[2], "posterior samples returned \n") else cat("Transmission tree inferred, but posterior samples not returned\n") 
}

# FUNCTION summary.epinet
# Summary method for class epinet

summary.epinet <- function(object, ...)
{
	cat("Epidemic parameter posterior summaries: \n")
	print(summary(data.frame("Beta" = object$beta, "Theta_e" = object$thetae, "Theta_i" = object$thetai, "k_e" = object$ke, "k_i" = object$ki)))
	cat("Network parameter posterior summaries: \n")
	print(summary(object$eta))
}

# FUNCTION plot.epinet
# Plot method for class epinet

plot.epinet <- function(x, index = dim(x$transtree)[2], lwd = 1,
leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, main = "Transmission Tree",
xlab = "Time", ylab= "", e.col = "black", i.col = "red", lty.transmission = 3,
marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...) {
    
    if(is.null(x$transtree) || is.null(index)) stop("Error: No inferred transmission tree samples found to plot.")
    
    epi <- buildepifromoutput(x, index)
    
    plotepitree(epi, lwd = lwd, leaf.labs = leaf.labs, leaf.cex = leaf.cex,
    zero.at.start = zero.at.start, main = main, xlab = xlab, ylab= ylab, e.col = e.col,
    i.col = i.col, lty.transmission = lty.transmission, marktransitions = marktransitions,
    label.trans = label.trans, cex.trans = cex.trans, ...)
}
