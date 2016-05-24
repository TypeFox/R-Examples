phtMCMC = function(x, states, beta, nu, zeta, n, mhit=1, resume=NULL, silent=FALSE) {
	# Fake what we need for code reuseability
	method <- "MHRS"
	censored <- rep(FALSE, length(x))
	
	# Note, nu is assumed to specify rowwise the prior parameters
	if(length(nu) != states*states)
		stop("nu must specify one prior Gamma shape parameter per element of the Phase-type generator matrix\n    [Hint: so to fit a Phase-type with m non-absorbing states, nu must have m*m elements]")
	if(length(zeta) != states)
		stop("zeta must specify one prior Gamma reciprocal scale parameter per non-absorbing row of the Phase-type generator matrix\n    [Hint: so to fit a Phase-type with m non-absorbing states, zeta must have m elements]")
	
	# Check n and mhit
	if(!is.numeric(n) || as.integer(n)<1) stop(n, " is an invalid number of MCMC iterations.\n    [Hint: must be >= 1]\n\n")
	if(!is.numeric(mhit) || as.integer(mhit)<0) stop(mhit, " is an invalid number of Metropolis-Hastings iterations.\n    [Hint: must be >= 1]\n\n")
	
	# Form T
	T <- matrix(paste("S", 1:states, rep(1:states, each=states), sep = ""), nrow=states, ncol=states)
	T <- cbind(T, paste("s", 1:states, sep = ""))
	T <- rbind(T, rep(0, states+1))
	diag(T) <- 0
  
	# Check beta
	if(length(beta) != states) stop("beta should be a vector of length ", states, " for the generator specified.")
	if(sum(beta<0) > 0) stop("beta is not a valid parameter of a Dirichlet distribution.\n    [Hint: every element must be non-negative]\n\n")
	
	# Create names for nu and zeta (after making zeta the right size to reuse C code)
	nu <- as.list(nu)
	names(nu) <- setdiff(c(t(T)), c("0"))
	zeta <- as.list(rep(zeta, each=states+1))
	names(zeta) <- setdiff(c(t(T)), c("0"))
	
	# Ensure ordered by var name
	nu <- nu[sort(names(nu))]
	zeta <- zeta[sort(names(zeta))]
	
	# Extract variable names
	varNames <- sort(unique(c(T)))
	varNames <- subset(varNames, varNames!="0")
	
	# Check names of nu params match
	nuNames <- sort(names(nu))
	if(length(nuNames)!=length(varNames) || sum(nuNames==varNames)-length(varNames)!=0)
		stop("variables specified in matrix don't match those in prior nu\n    [Hint: there should be one nu entry for every different variable being estimated ... so nu should be a list in R with each value's tag being the variable name matching that in the matrix of variables]\n\n")
	
	# Check names of zeta params match
	zetaNames <- sort(names(zeta))
	if(length(zetaNames)!=length(varNames) || sum(zetaNames==varNames)-length(varNames)!=0)
		stop("variables specified in matrix don't match those in prior zeta\n    [Hint: there should be one zeta entry for every different variable being estimated ... so zeta should be a list in R with each value's tag being the variable name matching that in the matrix of variables]\n\n")
	
	# Check that resume is well formed
	start <- -1
	if(!is.null(resume)) {
		if(!is.mcmc(resume)) stop("resume must be an mcmc object\n    [Hint: mcmc objects are defined in the coda package]\n\n")
		if(length(colnames(resume)) != length(varNames)) stop("the variables in resume do not match the variables in generator\n    [Hint: are you certain you are trying to resume the correct chain?]\n\n")
		for(i in 1:dim(resume)[2]) {
			if(colnames(resume)[i] != varNames[i])
				stop("the variable names in resume do not match the variable names in generator\n    [Hint: are you certain you are trying to resume the correct chain?]\n\n")
		}
		start <- resume[dim(resume)[1],]
		n <- n+1 # because otherwise we repeat resume value
	}
	
	# Create return data frame
	ret <- matrix(0, nrow=n, ncol=length(varNames))
	
	# Transform variable names to variable indicies for C
	dimT <- dim(T)[1]
	TN <- matrix(0, dimT, dimT)
	for(i in 1:dimT) {
		for(j in 1:dimT) {
			TN[i,j] <- which(c("0",varNames)==T[i,j])-1
		}
	}
	
	# Organize correct method identification in C
	methodKey <- data.frame(MHRS=1, ECS=2, DCS=4)
	if(length(setdiff(method, names(methodKey)))>0) {
		cat("Error: unknown sampling methods (", setdiff(method, names(methodKey)), ")\n    [Hint: Options are MHRS, ECS and DCS]\n\n")
		return()
	}
	methodNum <- sum(methodKey[unique(method)])
	
	# Run MCMC
	res <- .C("LJMA_Gibbs", it=as.integer(n), mhit=as.integer(mhit), method=as.integer(methodNum), n=as.integer(dimT-1), m=as.integer(length(varNames)), nu=as.double(nu), zeta=as.double(zeta), T=as.integer(TN), C=as.double(matrix(1.0, nrow=dimT, ncol=dimT)), y=as.double(x), l=as.integer(length(x)), censored=as.integer(censored), start=as.double(start), silent=as.integer(silent), res=as.double(ret))
	
	if(!is.null(resume)) {
		ret <- rbind(resume[-dim(resume)[1],], matrix(res$res, nrow=n, ncol=length(varNames)))
	} else {
		ret <- matrix(res$res, nrow=n, ncol=length(varNames))
	}
	ret <- as.data.frame(ret)
	names(ret) <- varNames
	
	#beta, nu, zeta, n, censored=rep(FALSE, length(x)), method=c("ForwardUnconditional", "ReverseUnconditional", "ForwardConditional", "ReverseConditional", "Hobolth"), mhit=1, resume=NULL, stoprule=NULL)
	res <- list(samples=as.mcmc(ret), data=x, vars=varNames, T=T, beta=beta, nu=nu, zeta=zeta, iterations=n, censored=censored, method=method, MHit=mhit)
	class(res) <- "phtMCMC"
	res
}
