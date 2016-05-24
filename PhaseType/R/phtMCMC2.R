phtMCMC2 = function(x, T, beta, nu, zeta, n, censored=rep(FALSE, length(x)), C=matrix(1.0, nrow=dim(T)[1], ncol=dim(T)[2]), method="ECS", mhit=1, resume=NULL, silent=FALSE) {
	# Ensure ordered by var name
	nu <- nu[sort(names(nu))]
	zeta <- zeta[sort(names(zeta))]
	
	# Check n and mhit
	if(!is.numeric(n) || as.integer(n)<1) stop(n, " is an invalid number of MCMC iterations.\n    [Hint: must be >= 1]\n\n")
	if(!is.numeric(mhit) || as.integer(mhit)<0) stop(mhit, " is an invalid number of Metropolis-Hastings iterations.\n    [Hint: must be >= 1]\n\n")
	
	# Check T well formed
	dimT <- dim(T)
	if(dimT[1]!=dimT[2]) stop("matrix of variables must be square\n    [Hint: every continuous-time Markov chain generator must be square]\n\n")
	dimT <- dimT[1]
	if(length(unique(diag(T))) != 1 || unique(diag(T)) != "0")
		stop("diagonal of matrix of variables must be zeros\n    [Hint: this is just because the diagonal of a continuous-time Markov chain generator is fully specified by the off-diagonal entries and therefore is not part of the inferential procedure]\n\n")
	if(length(unique(T[dimT,])) != 1 || unique(T[dimT,]) != "0")
		stop("last row of matrix of variables must represent absorbing state (and so be all zeros)\n    [Hint: this is just a formatting requirement ... the absorbing state can always be made the last state without loss of generality]\n\n")
	
	# Check beta
	if(length(beta) != dimT-1) stop("beta should be a vector of length ", dimT-1, " for the generator specified.")
  if(sum(beta<0) > 0) stop("beta is not a valid parameter of a Dirichlet distribution.\n    [Hint: every element must be non-negative]\n\n")
	
	# Check C well formed
	if(sum(dim(C)==dimT) < 2)
		stop("dimension of C must match dimension of T\n	[Hint: C is the matrix of constant multiples of parameter values]\n\n")
	
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
	res <- .C("LJMA_Gibbs", it=as.integer(n), mhit=as.integer(mhit), method=as.integer(methodNum), n=as.integer(dimT-1), m=as.integer(length(varNames)), nu=as.double(nu), zeta=as.double(zeta), T=as.integer(TN), C=as.double(C), y=as.double(x), l=as.integer(length(x)), censored=as.integer(censored), start=as.double(start), silent=as.integer(silent), res=as.double(ret))
	
	if(!is.null(resume)) {
		ret <- rbind(resume[-dim(resume)[1],], matrix(res$res, nrow=n, ncol=length(varNames)))
	} else {
		ret <- matrix(res$res, nrow=n, ncol=length(varNames))
	}
	ret <- as.data.frame(ret)
	names(ret) <- varNames
	
	res <- list(samples=as.mcmc(ret), data=x, vars=varNames, T=T, beta=beta, nu=nu, zeta=zeta, iterations=n, censored=censored, method=method, MHit=mhit)
	class(res) <- "phtMCMC"
	res
}
