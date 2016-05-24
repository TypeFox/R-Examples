
# File NetworkEpiBayesSEIR.R

# FUNCTION itimestartvalues
# Internal function used by epibayesmcmc
# Determines start values for Infective times (I) if needed

itimestartvalues <- function(epidata)
{	
	# Want to scroll through in order of increasing E times -- indexset gives this order
	indexset <- order(rank(epidata[,3]))
	
	# Set initial I time
	epidata[indexset[1],4] <- runif(1,epidata[indexset[1],3],epidata[indexset[2],3])
	currrectime <- nextrectime <- epidata[indexset[1],5]
	nextrec <- indexset[1]
	
	# Cycle through the remaining individuals
	for (i in 2:length(epidata[,4]))
	{		
		if (epidata[indexset[i],3] > currrectime)
		{
			currrectime <- nextrectime
			epidata[nextrec,4] <- runif(1,epidata[nextrec,3],epidata[indexset[i],3])
		}		
		epidata[indexset[i],4] <- runif(1,epidata[indexset[i],3],epidata[indexset[i],5])
		if(epidata[indexset[i],5] > nextrectime)
		{
			nextrectime <- epidata[indexset[i],5]
			nextrec <- indexset[i]
		}
	}
	
	return(epidata[,4])
}


# FUNCTION etimestartvalues
# Internal function used by epibayesmcmc
# Determines start values for Exposure times (E) if needed

etimestartvalues <- function(epidata,initialoffset)
{	
	for (i in 1:length(epidata[,3]))
	{
		currpar <- -999
		maxrand <- 0
		posspar <- 0
		for (j in 1:length(epidata[,3]))
		{
			if (epidata[i,4] > epidata[j,4])	
			posspar <- runif(1,0,1)
			if (posspar > maxrand)
			{
				currpar <- j
				maxrand <- posspar
			}
		}
		if (currpar == -999) 
			epidata[i,3] <- epidata[i,4] - initialoffset
		else
			epidata[i,3] <- runif(1,epidata[currpar,4], min(epidata[currpar,5],epidata[i,4]))
	}	
	return(epidata[,3])
}


# FUNCTION epibayesmcmc
# R wrapper for C function that performs Bayesian MCMC inference

epibayesmcmc <- function(epidata, dyadiccovmat, nsamp, thinning, bprior, tiprior, 
    teprior, etaprior, kiprior, keprior, etapropsd, priordists = "gamma", betapriordist = priordists, 
    thetaipriordist = priordists, thetaepriordist = priordists, etapriordist = rep("normal", times=etapars), 
    kipriordist = priordists, kepriordist = priordists, extrathinning = FALSE, inferEtimes = FALSE, 
    inferItimes = FALSE, parentprobmult = 1, verbose = TRUE, burnin = 0)
{		
	
	# Check overall input format
	
	if (!is.numeric(epidata)) stop("Invalid input format: epidata must by numeric")
	if (length(dim(epidata)) != 2) stop ("Invalid input format: epidata must be a two-dimensional matrix")
		
	# N = number of individuals in population
	
	N = dim(epidata)[1] 

	# Do some processing on dyadic covariate matrix

	etapars <- dim(dyadiccovmat)[2] - 2
	
	if (length(etapropsd) != etapars) stop("Invalid input: length(etapropsd) must equal the number of eta parameters in the model")
		
	if( is.null(dimnames(dyadiccovmat)[[2]]))
		etanames <- paste("factor",1:etapars)
	else etanames <- dimnames(dyadiccovmat)[[2]][-(1:2)]

	# Do some cleanup of parameters
	
	nsamp <- floor(nsamp)
	thinning <- floor(thinning)
	extrathinning <- floor(extrathinning)
	
	maxmove <- 11
		
	# Do some processing on infection/recovery times:
	
	# (1) Shift infection and recovery times so that the first recovery happens at time 0
	# 	 This is necessary since the prior on the initial infected is from -Inf to 0
	
	epidata[,3:5] <- epidata[,3:5] - min(epidata[,5], na.rm = TRUE)
			
	# (2a) Set initial values of parameters
	# 	For beta, thetai, thetae, ki, ke, we are starting these parameters at the mean of their respective prior distributions
	#		(except for the theta parameters -- we are starting them at their modes, since they may not have means)
	# 	Beta, ki, ke, thetai, and thetae have Gamma/IG (or uniform) priors
	#	Starting all of the eta parameters at 0 (for flat priors) or the prior mean (for normal priors)
	
	if(betapriordist == "gamma") initbeta <- bprior[1]*bprior[2] else initbeta <- (bprior[1]+bprior[2])/2
	if(thetaipriordist == "gamma") initthetai <- tiprior[2]/(tiprior[1] + 1) else initthetai <- (tiprior[1]+tiprior[2])/2
	if(thetaepriordist == "gamma") initthetae <- teprior[2]/(teprior[1] + 1) else initthetae <- (teprior[1]+teprior[2])/2
	if(kipriordist == "gamma") initki <- kiprior[1]*kiprior[2] else initki <- (kiprior[1]+kiprior[2])/2
	if(kepriordist == "gamma") initke <- keprior[1]*keprior[2] else initke <- (keprior[1]+keprior[2])/2

	initeta <- rep(0,times=etapars)
	for (i in 1:etapars)
		if(etapriordist[i] == "normal")
			initeta[i] <- etaprior[2*i-1]
				
	if (verbose)
	{
		cat("epinet run started at ",format(Sys.time(), "%Y-%m-%d %X"),"\n")
		cat("\n Initial parameter values \n")
		cat("------------------------------ \n")
		cat("Beta = ",initbeta,"\n")
		for (i in 1:etapars)
			cat("eta_",i-1, " = ", initeta[i], "\n")
		cat("Theta_I = ",initthetai,"\n")
		cat("Theta_E = ",initthetae,"\n")
		cat("k_I = ",initki,"\n")
		cat("k_E = ",initke,"\n \n")
	}
		
	# (2b) Set values for prior distributions
	#	0 indicates a uniform prior dist
	#	1 indicates a gamma/IG prior dist
	#	2 indicates a flat prior (not currently used)
	#	3 indicates a normal prior (mean, sd)

	bpriordistnum <- 1*(betapriordist == "gamma")
	tipriordistnum <- 1*(thetaipriordist == "gamma")
	tepriordistnum <- 1*(thetaepriordist == "gamma")
	kipriordistnum <- 1*(kipriordist == "gamma")
	kepriordistnum <- 1*(kepriordist == "gamma")	
	etapriordistnum <- 2 + 1*(etapriordist == "normal")
	
	# (3a) Find number of infecteds in population
	
	ninf <- min(which(is.na(epidata[,5])),dim(epidata)[1] + 1) - 1
	
	# Set E, I, and R times of susceptibles to "infinity"; basically, any constant that is larger than all of the other times (by 10 here)
	if (N > ninf) epidata[(ninf+1):(N),3:5] <- max(epidata[1:ninf,(3:5)],na.rm=TRUE) + 10
	
	if (verbose) cat("Epidemic data: ", ninf, "infecteds,", N - ninf, "susceptibles, ", N, "total individuals. \n \n")
	
	# (3b) Determine starting values for E/I times, if we are inferring these times (if not, make sure we have legitimate values)
	
	if (inferItimes) 
	{
		if (inferEtimes) 
			epidata[(1:ninf),4] <- epidata[(1:ninf),5] - rgamma(ninf,initki,1/initthetai) 
		else 
		{
			if (sum(is.na(epidata[(1:ninf),3])) > 0) stop("Invalid input data: NAs found in Exposure times -- run with InferEtimes = TRUE.") 
			if (sum(epidata[(1:ninf),3] > epidata[(1:ninf),5]) > 0) stop("Invalid input data: Exposure times cannot occur after Removal times.")
			epidata[(1:ninf),4] <- itimestartvalues(epidata[(1:ninf),]) 
		}
	}
	else
		if (sum(is.na(epidata[(1:ninf),4])) > 0) stop("Invalid input data: NAs found in Infective times -- run with InferItimes = TRUE.") 
		
	if (inferEtimes) 
	{
		if (sum(epidata[(1:ninf),4] > epidata[(1:ninf),5]) > 0) stop("Invalid input data: Infective times cannot occur after Removal times.")
		epidata[(1:ninf),3] <- etimestartvalues(epidata[(1:ninf),],initke*initthetae)	
	}
	else
		if (sum(is.na(epidata[(1:ninf),3])) > 0) stop("Invalid input data: NAs found in Exposure times -- run with InferEtimes = TRUE.") 
	
	# (4) Set up storage for extra stuff returned (E/I times and transmission tree)
	#	Note: if we are not inferring the E/I times, then we just return NULL for these values
	#	There are generally too many times to return each iteration (we run out of space)
	#	So instead, if we want to return these, we use the variable "extrathinning" as the extra thinning interval

    numsamp <- floor(nsamp/thinning)
    if (extrathinning > 0) numsamptimes <- floor(nsamp/(thinning*extrathinning)) else numsamptimes <- NULL
	#cat(numsamp, " ", numsamptimes, "\n")

	if (extrathinning == 0)
	{
		storeexptimes <- NULL
		storeinftimes <- NULL
		storetranstree <- NULL
	} else
	{
		storeexptimes <- array(0,ninf*numsamptimes)
		storeinftimes <- array(0,ninf*numsamptimes)
		storetranstree <- array(0,ninf*numsamptimes)
	#cat(dim(storeexptimes))
	}
    
    
	# (5) Do translation (renaming) of nodes and probable parent info
	epidata[,1:2] <- floor(epidata[,1:2])
	
	# (5a) Check for invalid parent references
	for (i in 1:ninf)
		if (!is.na(epidata[i,2]) & (!(epidata[i,2] %in% epidata[(1:ninf),1]) | (epidata[i,2] == epidata[i,1]))) 
			stop("Invalid input data: bad parent reference")
	
	# (5b) Fill in any missing IDs	
	nextid <-max(epidata[,1],0,na.rm=TRUE)+1
	for (i in 1:length(epidata[,1]))
		if (is.na(epidata[i,1]))
		{
			epidata[i,1] <- nextid
			nextid <- nextid +1
		}
		
	# (5c) Check for duplicate IDs	
	if(length(epidata[,1]) > length(unique(epidata[,1])))
		stop("Invalid input data: non-unique node IDs")
		
	# (5d) Rename Node IDs
	parentprior <- rep(0, times=ninf)
	for (i in 1:ninf)
		if(!is.na(epidata[i,2]))
			parentprior[i] <- which(epidata[(1:ninf),1] == epidata[i,2])		
			
	# (5e) Re-sort dyadic covariates

	counter <- 1	
	sorteddyadiccovmat <- matrix(ncol=etapars,nrow=N*(N-1)/2)

	for (i in 1:(N-1))
	{
		for (j in (i+1):N)
		{
			if (epidata[i,1] < epidata[j,1])
			{
				minindex <- epidata[i,1]
				maxindex <- epidata[j,1]
			} else
			{
				minindex <- epidata[j,1]
				maxindex <- epidata[i,1]
			}			
			firstindex <- which(dyadiccovmat[,1]==minindex)
			secondindex <- which(dyadiccovmat[,2]==maxindex)
			sorteddyadiccovmat[counter,] <- dyadiccovmat[intersect(firstindex,secondindex),-(1:2)]
			counter <- counter + 1
		}
	}
	
	changed <- FALSE
	for (i in 1:etapars)
		for(j in 1:(N*(N-1)/2))
			if (dyadiccovmat[j,i+2] != sorteddyadiccovmat[j,i])
				changed <- TRUE

	if (verbose)
	{
		if (changed) cat("Re-sorted dyadic covariates.\n")
			else cat("Dyadic covariates not re-sorted.\n")
		cat("Coefficient names: \n")
		for (i in 1:etapars)
			cat("eta", i-1, ": ", etanames[i], "\n")
		cat("\n")
	}

	# Call C function to do actual MCMC routine
    
	output <- .C("epigraphmcmcc",as.double(epidata[,3]),as.double(epidata[,4]), as.double(epidata[,5]), as.integer(etapars),
		as.double(sorteddyadiccovmat), as.integer(nsamp), as.integer(thinning),as.double(bprior), as.double(tiprior), as.double(teprior), as.double(etaprior), 
		as.double(kiprior), as.double(keprior), as.integer(ninf), as.integer(N), as.double(initbeta),as.double(initthetai), as.double(initki), as.double(initthetae), 
		as.double(initke), as.double(initeta), as.integer(bpriordistnum), as.integer(tipriordistnum), as.integer(tepriordistnum), as.integer(kipriordistnum), 
		as.integer(kepriordistnum) ,as.integer(etapriordistnum), as.double(etapropsd), accept=as.integer(array(0,maxmove-3)), propose=as.integer(array(0,maxmove-3)), 
		llkd = as.double(array(0,numsamp)), betaout = as.double(array(0,numsamp)),
		thetaiout=as.double(array(0,numsamp)), kiout=as.double(array(0,numsamp)), thetaeout=as.double(array(0,numsamp)),
		keout=as.double(array(0,numsamp)), eta=as.double(array(0,numsamp*etapars)), initexp=as.integer(array(0,numsamp)),
		initexptime=as.double(array(0,numsamp)), exptimes = as.double(storeexptimes), inftimes = as.double(storeinftimes),
		transtree = as.integer(storetranstree), as.integer(extrathinning), as.integer(1*inferEtimes), as.integer(1*inferItimes), as.integer(parentprior), 
		as.integer(parentprobmult), as.integer(1*verbose), as.integer(burnin), as.integer(numsamp), as.integer(numsamptimes)
        #	, PACKAGE="epinet"
	)
	
	# Post-processing	
	
	# Arrange eta into an array, with the samples for each group corresponding to one column
	etaout <- aperm(array(output$eta,dim=c(etapars,numsamp)))
	dimnames(etaout) <- list(NULL,covariates=etanames)
	
	# Reverse translate nodes for initinf (back to original node IDs)
	initexpout <- epidata[output$initexp,1]
	
	if(extrathinning == 0)
	{	
		expout <- NULL
		infout <- NULL
		transtreeout <- NULL
	} else
	{
		expout <- array(output$exptimes,dim=c(ninf,numsamptimes))
		infout <- array(output$inftimes,dim=c(ninf,numsamptimes))
		# Back translate transmission tree to their original node IDs
		output$transtree[output$transtree == -999] <- NA
		transtreeout <- array(epidata[output$transtree,1],dim=c(ninf,numsamptimes))
	}		
	
	return(list(accept=output$accept,propose=output$propose,llkd=output$llkd, beta=output$betaout,thetai=output$thetaiout,thetae=output$thetaeout,ki=output$kiout,
		ke=output$keout, eta=etaout,initexp=initexpout, initexptime=output$initexptime,exptimes=expout,inftimes=infout,rectimes=epidata[,5],nodeid=epidata[,1],
		transtree=transtreeout))
	
}
