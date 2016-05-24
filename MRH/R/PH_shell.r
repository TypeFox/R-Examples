###########################################################################################
# This code is used to estimate the parameters of a given data set using the MRH
#	methodology.  One to many covariates can be included under the proportional 
#	hazards assumption.
###########################################################################################

PHshell = function(Mval, Ti, delta, X, outfilename, censortime,
prune.indc, burnIn, maxIter, thin, RmpInit, a.init, lambda.init, beta.init, 
k, gamma.mp, GR, fix.burnIn, fix.thin, fix.max,
writenum, sysseconds, systemtime, checknum, sampletype, continue.chain, k.init,
gamma.init, loopctr.start){

    ###########################################################################
	#		If X matrix is available, standardize the continuous 
	#		variables and label the betas
	###########################################################################
	if(!is.null(X)){
		n = length(Ti)
		betaNames = colnames(X)
		numParams = ncol(X)
				
		# Find the bounds on the beta estimates, used for sampling
		betaInfo = summary(coxph(Surv(Ti, delta) ~ as.matrix(X)))$coeff
		if(length(which(betaInfo[,1] < -50 | betaInfo[,1] > 50)) > 0){
			stop(paste("Algorithm will not converge.  Covariate estimate for variable(s)", 
					   which(betaInfo[,1] < -50 | betaInfo[,1] > 50), "is out of range."))
		}
		betaLB = betaInfo[,1]-betaInfo[,3]*5*sqrt(n)
		betaUB = betaInfo[,1]+betaInfo[,3]*5*sqrt(n)
		betaLB[betaLB < -50] = -50
		betaUB[betaUB > 50] = 50
		betaInfo[,1][betaInfo[,1] < betaLB | betaInfo[,1] > betaUB] = 0
			
		# Check the bounds on the continuous covariates 
		standardizeinfo = prepX.contin(X, betaLB, betaUB)
		Xstdz = as.matrix(standardizeinfo$X)
		stdIndex = standardizeinfo$stdIndex
		Xmeans = rep(0, numParams)
		Xmeans[stdIndex] = standardizeinfo$Xmeans
		Xsds = rep(1, numParams)
		Xsds[stdIndex] = standardizeinfo$Xsds
		X = as.matrix(X)
		
        betas.LB.stdz = betas.UB.stdz = NULL
		if(length(stdIndex) > 0){
			betas.LB.stdz = betaLB*Xsds
			betas.UB.stdz = betaUB*Xsds
			betas.LB.stdz[betas.LB.stdz < -50] = -50
			betas.UB.stdz[betas.UB.stdz > 50] = 50
		}
	} else {
		numParams = 0
		n = length(Ti)
		Xmeans = Xsds = stdIndex = NULL
		X = matrix(1, ncol = 1, nrow = n)
		beta = matrix(0, ncol = 1, nrow = 1)
	}
    mcmcInfoMatrix = as.data.frame(matrix(c(censortime, n, outfilename), nrow = 1))
    names(mcmcInfoMatrix) = c('maxStudyTime', 'n', 'outfolder')
    write.table(mcmcInfoMatrix, paste(outfilename, '/MCMCInfo.txt', sep = ''), row.names = FALSE)

	#########################################################
	# Create all variables that only need to be set up once 
	# during the estimation.
	#########################################################
	# mat01 is the PI matrix on p328 of the SINICA paper
	mat01 = calc_mat01(M = Mval)
	# TiBR holds the bin number for each persons time, as well
	# as the ratio of time each person spent in the final bin.	
	TiBR = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)

	# Calculate temp, which holds the indicators for both the Rmp and 1-Rmp values.
	# These indicators are needed to calculate the Rmp posteriors.
	temp = calcXIndics(mat01 = mat01, bin_ratio = TiBR, M = Mval)
	# RmpInd holds the indicator for the Rmp values, with the columns as follows: 
	# {R10, R20, R21, R30, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	RmpInd = matrix(temp[,1:(2^Mval-1)], ncol = 2^Mval-1)
	# one_RmpInd holds the 1-Rmp indicators, with the columns as follows:
	# {1-R10, 1-R20, 1-R21, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	one_RmpInd = matrix(temp[,-c(1:(2^Mval-1))], ncol = 2^Mval-1)
	
	# The inBin and failBin matrices are needed for calcluating H and F quickly.  
	whereFail = TiBR$bin+1
	whereFail[whereFail > 2^Mval] = 2^Mval
	# inBin is a matrix that holds a '1' for every bin the subject 
	# survived through, and the ratio in the bin they failed in, with
	# zeroes in all columns beyond the failure bin.
	# failBin is a matrix that holds a '1' in the column 
	# corresponding to the bin that the subject failed in, and '0' otherwise.
	inBin = failBin = matrix(0, ncol = 2^Mval, nrow = n)
	inBin[,1] = TiBR$ratio
	failBin[whereFail == 1, 1] = 1
	for(j in 2:2^Mval){
		inBin[whereFail == j, 1:(j-1)] = 1
		inBin[whereFail == j, j] = TiBR$ratio[whereFail == j]
		failBin[whereFail == j, j] = 1
	}
	
	# mvec will hold the value of m for m = 1,...,M, with
	# each number repeated 2^(m-1) times.  This is used in 
	# the calculations of the Rmp posteriors.
	mvec = 1
	if(Mval > 1){ for(mctr in 2:Mval){mvec = c(mvec, rep(mctr, 2^(mctr-1)))}}
		
	#########################################################
	# Initialize parameters
	#########################################################
	Rmp = rep(0.5, 2^Mval-1)
	a = a.init
	lambda = -log(mean(delta))/a
	if(!is.null(lambda.init)){	lambda = lambda.init	}
	
	# If user has specified their own Rmp initialization values, then use them
	if(!is.null(RmpInit)){ Rmp = RmpInit	}
	if(numParams > 0){
		if(is.null(beta.init)){	beta = betaInfo[,1]
		} else {	beta = beta.init	}
	} 
	# If the user is not using pruned bins, set Rmp sampling index vector so that 
	# all Rmps are sampled.  If the user would like to use pruned bins, set the 
	# Rmp sampling index vector so that those Rmps are not sampled or even checked.
	RmpNoPruneIndx = 1:(2^Mval-1)
	if(!is.null(prune.indc)){  RmpNoPruneIndx = RmpNoPruneIndx[-which(prune.indc == 1)]	}
	# If user is going to use the Gelman-Rubin test statistic, then jiggle the
	# initial values to cover the range for each parameter.
	if(GR == TRUE & continue.chain == FALSE){
		Rmp[RmpNoPruneIndx] = runif(length(RmpNoPruneIndx), 0, 1)
		a = round(runif(1, 0, 10))
		lambda = -log(mean(delta))/a
		if(numParams > 0){  for(i in 1:numParams){  beta[i] = runif(1, betaLB[i], betaUB[i])	}}
	}
	# Record all initialized values
    # Make the labels for the Rmps so that the output is easier to read
    RmpRowNums = NULL
    for(mctr in 1:Mval){ for(pctr in 0:(2^(mctr-1)-1)){
        RmpRowNums = c(RmpRowNums, paste(mctr, '.', pctr, sep = ''))
    }}
	tempRmp = as.data.frame(matrix(Rmp, nrow = 1))
	names(tempRmp) = paste('Rmp',RmpRowNums, sep = '')
	if(numParams > 0){
		tempbetas = as.data.frame(matrix(beta, nrow = 1))
		names(tempbetas) = betaNames	
	} else {	tempbetas = NULL	}
	initialValues = list(Rmp.init = tempRmp, beta.init = tempbetas, a.init = a, lambda.init = lambda)
	
    if(sampletype == 'simple'){
        output = PH_samplesimple(Mval, delta, X, Xstdz, outfilename, burnIn, maxIter, thin, k, gamma.mp,
                                    GR, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime,
                                    checknum, n, betaNames, numParams, betaLB, betaUB, beta,
                                    betas.LB.stdz, betas.UB.stdz, stdIndex, Xmeans, Xsds, mat01,
                                    RmpInd, one_RmpInd, inBin, failBin, mvec, Rmp, a, lambda, RmpNoPruneIndx,
                                    RmpRowNums, initialValues, continue.chain, loopctr.start)
    } else if(sampletype == 'konly'){
        output = PH_samplek(Mval, delta, X, Xstdz, outfilename, burnIn, maxIter, thin, gamma.mp,
                                    GR, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime,
                                    checknum, n, betaNames, numParams, betaLB, betaUB, beta,
                                    betas.LB.stdz, betas.UB.stdz, stdIndex, Xmeans, Xsds, mat01,
                                    RmpInd, one_RmpInd, inBin, failBin, mvec, Rmp, a, lambda, RmpNoPruneIndx,
                                    RmpRowNums, initialValues, continue.chain, k.init, loopctr.start)
    } else if(sampletype == 'gonly'){
        output = PH_samplegamma(Mval, delta, X, Xstdz, outfilename, burnIn, maxIter, thin, k,
                                    GR, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime,
                                    checknum, n, betaNames, numParams, betaLB, betaUB, beta,
                                    betas.LB.stdz, betas.UB.stdz, stdIndex, Xmeans, Xsds, mat01,
                                    RmpInd, one_RmpInd, inBin, failBin, mvec, Rmp, a, lambda, RmpNoPruneIndx,
                                    RmpRowNums, initialValues, continue.chain, gamma.init, loopctr.start)
    } else {
        output = PH_samplekandgamma(Mval, delta, X, Xstdz, outfilename, burnIn, maxIter, thin, 
                                    GR, fix.burnIn, fix.thin, fix.max, writenum, sysseconds, systemtime,
                                    checknum, n, betaNames, numParams, betaLB, betaUB, beta,
                                    betas.LB.stdz, betas.UB.stdz, stdIndex, Xmeans, Xsds, mat01,
                                    RmpInd, one_RmpInd, inBin, failBin, mvec, Rmp, a, lambda, RmpNoPruneIndx,
                                    RmpRowNums, initialValues, continue.chain, k.init, gamma.init, loopctr.start)
    }
   
	print("Estimation routine finished, preparing results....", quote = FALSE)
	
	outputdata = list(assessIndex = output$assessIndex, numberofEstParams = 2 + 1 + numParams + length(RmpNoPruneIndx),
                        numParams = numParams, failBin = failBin, inBin = inBin, GR = GR, loopctr = output$loopctr,
                        convergence = output$convergence, initialValues = output$initialValues, burnIn = output$burnIn, thin = output$thin)
    return(outputdata)
	
	
}
