###########################################################################################
# Created 6Sep13
###########################################################################################

NPHBAshell = function(Mval, Ti, Xfixed, XNPH, delta, outfilename,
prune.indc, burnIn, maxIter, thin, RmpInit, a.init, k, gamma.mp, 
lambda.init, beta.init, censortime, GR, convGraphs, fix.burnIn, fix.thin, fix.max, 
writenum, sysseconds, systemtime, checknum, sampletype, continue.chain, k.init,
gamma.init, loopctr.start){
    
	###########################################################################
	#	Standardize the continuous variables and label the betas
	###########################################################################
	n = length(Ti)
    mcmcInfoMatrix = as.data.frame(matrix(c(censortime, n, outfilename), nrow = 1))
    names(mcmcInfoMatrix) = c('maxStudyTime', 'n', 'outfolder')
    write.table(mcmcInfoMatrix, paste(outfilename, '/MCMCInfo.txt', sep = ''), row.names = FALSE)
			
	# Convert categorical variables to 0/1 
	namesHazGroups = names(table(XNPH))
	betaNPHnames_preBin = paste(colnames(XNPH), namesHazGroups[-1], sep = '.') 
			
	# Find the number of PH and NPH parameters, and create the beta names	
	if(!is.null(Xfixed)){	numPHParams = ncol(Xfixed)	} else {	numPHParams = 0	}
	numNPHParams = ncol(XNPH)
	numParams = numPHParams + numNPHParams
	if(is.na(betaNPHnames_preBin)[1]){ betaNPHnames_preBin = paste(rep('var', numNPHParams), 1:numNPHParams, sep = '') }
	betaNames = paste(rep(paste(betaNPHnames_preBin, '.bin', sep = ''), each = 2^Mval), 1:2^Mval, sep = '')
	if(numPHParams > 0){	betaNames = c(colnames(Xfixed), betaNames)		}
	# Find the bounds on the beta estimates, used for sampling
	binFail = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)$bin + 1
	binFail[binFail > 2^Mval] = 2^Mval
	binvar = factor(binFail, levels = 1:2^Mval)
	if(numPHParams > 0){
		betaInfo = summary(coxph(Surv(Ti, delta) ~ as.matrix(Xfixed)))$coeff
		if(length(which(betaInfo[,1] < -25 | betaInfo[,1] > 25)) > 0){
			stop(paste("Algorithm will not converge.  Covariate estimate for variable(s)", 
					   betaNames[which(betaInfo[,1] < -25 | betaInfo[,1] > 25)], "is (are) out of range."))
		}
		betaLB = betaInfo[,1]-betaInfo[,3]*3*sqrt(n)
		betaUB = betaInfo[,1]+betaInfo[,3]*3*sqrt(n)
		betaLB[betaLB < -25 | is.na(betaInfo[,4])] = -25
		betaLB[betaLB > 0] = -1
		betaUB[betaUB > 25 | is.na(betaInfo[,4])] = 25
		betaUB[betaUB < 0] = 1
		betaInfo[is.na(betaInfo[,4]),1] = 0
	}
	# Check the bounds on the continuous covariates 
	if(numPHParams > 0){
		stdIndex = NULL
		Xmeans = rep(0, numPHParams)
		Xsds = rep(1, numPHParams)
		standardizeinfo = prepX.contin(Xfixed, betaLB[1:numPHParams], betaUB[1:numPHParams])
		Xstdz = as.matrix(standardizeinfo$X)
		stdIndex = standardizeinfo$stdIndex
		Xmeans[stdIndex] = standardizeinfo$Xmeans
		Xsds[stdIndex] = standardizeinfo$Xsds
	}
			
	#########################################################
	# Create all variables that only need to be set up once 
	# during the estimation.
	#########################################################
	# mat01 is the PI matrix on p328 of the SINICA paper
	mat01 = calc_mat01(M = Mval)
	# Holds the sum of the columns of mat01 in pairs of two, and is used in the Rmpj posterior.
	newmat01 = cbind(1, mat01[,1:(ncol(mat01)/2-1)])
	# TiBR holds the bin number for each persons time, as well
	# as the ratio of time each person spent in the final bin.	
	TiBR = getbin.and.ratio(TJ = censortime, Ti = Ti, M = Mval)
			
	# The inBin and failBin matrices are needed for estimating betaj parameters.  
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
			
	# Calculate temp, which holds the indicators for both the Rmp and 1-Rmp values.
	# These indicators are needed to calculate the Rmp posteriors.
	temp = calcXIndics(mat01 = mat01, bin_ratio = TiBR, M = Mval)
	# RmpInd holds the indicator for the Rmp values, with the columns as follows: 
	# {R10, R20, R21, R30, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	RmpInd = matrix(temp[,1:(2^Mval-1)], ncol = 2^Mval-1)
	# one_RmpInd holds the 1-Rmp indicators, with the columns as follows:
	# {1-R10, 1-R20, 1-R21, etc}. The index goes from m = 1 to M, and p = 0 to 2^(m-1)-1.
	one_RmpInd = matrix(temp[,-c(1:(2^Mval-1))], ncol = 2^Mval-1)
	# mvec will hold the value of m for m = 1,...,M, with 
	# each number repeated 2^(m-1) times.  This is used in 
	# the calculations of the Rmp posteriors.
	mvec = 1
	if(Mval > 1){for(mctr in 2:Mval){	mvec = c(mvec, rep(mctr, 2^(mctr-1)))	}}
	# ind and one_ind are also needed in the calculation of the Rmp posteriors.
	ind = matrix(mat01[,1:(sum(2^(1:Mval))/2)*2-1], ncol = 2^Mval-1)
	one_ind = matrix(mat01[,1:(sum(2^(1:Mval))/2)*2], ncol = 2^Mval-1)
			
	# IDtypes will hold the types of different hazards #
	IDtypes = names(table(XNPH))
	IDtypes = IDtypes[order(IDtypes)]
	numHazards = length(IDtypes)
	indices = list(which(XNPH == IDtypes[1]))
	for(hazCtr in 2:numHazards){
		indices = c(indices, list(which(XNPH == IDtypes[hazCtr])))
	}
		
	#########################################################
	# Initialize parameters
	#########################################################
	Rmp = matrix(0.5, nrow = 2^Mval-1, ncol = numHazards)
	a = a.init
	if(length(a.init) < numHazards){	a = rep(a.init, numHazards)	}
	if(!is.null(lambda.init)){	
		if(length(lambda.init) < numHazards){	lambda = rep(lambda.init, numHazards)	
		} else {	lambda = lambda.init	}
	} else {	lambda = -log(sapply(1:numHazards, function(x) mean(delta[indices[[x]]]))[])/a	}
	# If user has specified their own Rmp initialization values, then use them
	if(!is.null(RmpInit)){ Rmp = RmpInit	}
	# Initialize the beta estimates
	if(is.null(beta.init) & !is.null(Xfixed)){	betas = betaInfo[,1]
	} else if(!is.null(beta.init)){	betas = beta.init	
	} else {	
		betas = 0	
		Xfixed = matrix(1, ncol = 1, nrow = n)
	}
	# If the user is not using pruned bins, set Rmp sampling index vector so that 
	# all Rmps are sampled.  If the user would like to use pruned bins, set the 
	# Rmp sampling index vector so that those Rmps are not sampled or even checked.
	RmpNoPruneIndx = NULL
	for(rmpctr in 1:numHazards){
		if(!is.null(prune.indc)){  
			RmpNoPruneIndx = c(RmpNoPruneIndx, list((1:(2^Mval-1))[which(prune.indc[,rmpctr] == 0)]))	
		} else {	RmpNoPruneIndx = c(RmpNoPruneIndx, list(1:(2^Mval-1)))	}
	}
	
	# If user is going to use the Gelman-Rubin test statistic, then jiggle the
	# initial values to cover the range for each parameter.
	if(GR == TRUE & continue.chain == FALSE){
		for(rmpctr in 1:numHazards){
			Rmp[RmpNoPruneIndx[[rmpctr]],rmpctr] = runif(length(RmpNoPruneIndx[[rmpctr]]), 0, 1)
		}
		a = runif(numHazards, 0, 10)
		lambda = -log(sapply(1:numHazards, function(x) mean(delta[indices[[x]]]))[])/a
        if(numPHParams > 0){
            for(i in 1:length(betas)){	betas[i] = runif(1, betaLB[i], betaUB[i])	}
        }
	}
			
	#########################################################
	# Start output file for parameter estimates
	#########################################################
	# Make the labels for the Rmps so that the output is easier to read	
	RmpRowNums = NULL
	for(mctr in 1:Mval){ for(pctr in 0:(2^(mctr-1)-1)){ 
		RmpRowNums = c(RmpRowNums, paste(mctr, '.', pctr, sep = '')) 
	}}
	RmpRowNums = paste(RmpRowNums, rep(namesHazGroups, each = 2^Mval-1), sep = '_')
	dnames = paste(paste('d', 1:2^Mval, sep = ''), rep(namesHazGroups, each = 2^Mval), sep = '_')

	# Record all initialized values
	tempRmpj = as.data.frame(matrix(Rmp, nrow = 1))
	names(tempRmpj) = paste('Rmp',RmpRowNums, sep = '')
	if(!is.null(Xfixed)){	
		tempbetas = as.data.frame(matrix(betas, nrow = 1))	
		names(tempbetas) = betaNames[1:numPHParams]
	} else {	tempbetas = NULL	}
	initialValues = list(Rmp.init = tempRmpj, beta.init = tempbetas, a.init = a, lambda.init = lambda)
    
    if(sampletype == 'simple'){
        output = NPHBA_simple(Mval, Xfixed, delta, outfilename, prune.indc, burnIn, maxIter, thin, Rmp, a, k,
                                gamma.mp, lambda, betas, fix.burnIn, fix.thin, fix.max, writenum, sysseconds,
                                systemtime, checknum, n, dnames, namesHazGroups, numPHParams, numParams,
                                betaNames, betaLB, betaUB, Xmeans, Xsds, Xstdz, stdIndex,
                                mat01, inBin, failBin, RmpInd, one_RmpInd, mvec, numHazards,
                                RmpNoPruneIndx, GR, initialValues, RmpRowNums, indices, continue.chain,
                                loopctr.start)
    } else if(sampletype == 'konly'){
        output = NPHBA_samplek(Mval, Xfixed, delta, outfilename, prune.indc, burnIn, maxIter, thin, Rmp, a,
                                gamma.mp, lambda, betas, fix.burnIn, fix.thin, fix.max, writenum, sysseconds,
                                systemtime, checknum, n, dnames, namesHazGroups, numPHParams, numParams,
                                betaNames, betaLB, betaUB, Xmeans, Xsds, Xstdz, stdIndex,
                                mat01, inBin, failBin, RmpInd, one_RmpInd, mvec, numHazards,
                                RmpNoPruneIndx, GR, initialValues, RmpRowNums, indices, continue.chain,
                                k.init, loopctr.start)
    } else if(sampletype == 'gonly'){
        output = NPHBA_samplegamma(Mval, Xfixed, delta, outfilename, prune.indc, burnIn, maxIter,
                                    thin, Rmp, a, k, lambda, betas, fix.burnIn, fix.thin, fix.max, writenum, sysseconds,
                                    systemtime, checknum, n, dnames, namesHazGroups, numPHParams, numParams,
                                    betaNames, betaLB, betaUB, Xmeans, Xsds, Xstdz, stdIndex,
                                    mat01, inBin, failBin, RmpInd, one_RmpInd, mvec, numHazards,
                                    RmpNoPruneIndx, GR, initialValues, RmpRowNums, indices, continue.chain,
                                    gamma.init, loopctr.start)
    } else {
        output = NPHBA_samplekandgamma(Mval, Xfixed, delta, outfilename, prune.indc, burnIn, maxIter,
                                        thin, Rmp, a, lambda, betas, fix.burnIn, fix.thin, fix.max, writenum, sysseconds,
                                        systemtime, checknum, n, dnames, namesHazGroups, numPHParams, numParams,
                                        betaNames, betaLB, betaUB, Xmeans, Xsds, Xstdz, stdIndex,
                                        mat01, inBin, failBin, RmpInd, one_RmpInd, mvec, numHazards,
                                        RmpNoPruneIndx, GR, initialValues, RmpRowNums, indices,
                                        continue.chain, k.init, gamma.init, loopctr.start)
    }
    
    # Number of estimated parameters: a, lambda, H for each hazard, number of unpruned Rmps, number of betas
	outputdata = list(assessIndex = output$assessIndex,
                        numberofEstParams = 3*numHazards + length(output$RmpAnalyzeIndex) + numPHParams,
                        numParams = numParams, failBin = failBin, inBin = inBin, GR = GR, loopctr = output$loopctr,
                        convergence = output$convergence, initialValues = output$initialValues, numPHParams = numPHParams,
                        namesHazGroups = namesHazGroups, burnIn = output$burnIn, thin = output$thin)

    return(outputdata)
	
}

