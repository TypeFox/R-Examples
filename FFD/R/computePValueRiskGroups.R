## Ian Kopacka
## 2011-07-04
##
## Function: computePValueRiskGroups
## 
## Computes the probability of finding no testpositives in
## a sample of a finite population where an imperfect test is used and the 
## population is stratified into risk groups. For each of the risk groups
## the population size, the sample size and the relative infection risk
## must be specified.
##
## Input parameters:
##     nPopulationVec...Integer vector. Population sizes of the risk groups.
##     nSampleVec.......Integer vector. Sample sizes of the risk groups.
##     nRelRiskVec......Numeric vector. (Relative) infection risks of the 
##                      risk groups.
##     nDiseased........Integer. Number of diseased in the total population 
##                      (according to the null hypothesis).
##     sensitivity......Numeric between 0 and 1. Sensitivity of test (diagnostic 
##                      test for one stage sampling, herd test for two stage 
##                      sampling).
##     specificity......Numeric between 0 and 1. Specificity of test (diagnostic 
##                      test for one stage sampling, herd test for two stage 
##                      sampling).
##
## Sources: 
## A.R. Cameron, F.C. Baldock, "A new probability fomula for 
## surveys to substantiate freedom from disease", Prev. Vet. Med. 34
## (1998), pp. 1 - 17.
## P.A.J.Martin, A.R. Cameron, M. Greiner, "Demonstrating freedom from disease
## using multiple complex data sources. !: A new methodology based on scenario
## trees", Prev. Vet. Med. 79 (2007), pp. 71 - 97.
##
## Calls:
##    computePValue.R
##    computeDiseaseMatrix.R
##
## Is called by:
##    computeOptimalSampleSize.R
##    computeAlphaLimitedSampling.R
##
computePValueRiskGroups <- function(nPopulationVec, nSampleVec, nRelRiskVec,
	nDiseased, sensitivity, specificity = 1){
	## Factor for normalization:
	theta <- sum(nPopulationVec*nRelRiskVec)
	## Vector of probabilities p_i:
	pVec <- nRelRiskVec*nPopulationVec/theta
	## Sensitivites:
	if (length(sensitivity) == 1) sensitivity <- rep(sensitivity, length(nPopulationVec))
	

	## Two Risk groups:
	if (length(nPopulationVec) == 2){
		## Possible number of diseased in RG 1:
		d1Vec <- max(c(0,nDiseased-nPopulationVec[2])):
			min(c(nDiseased,nPopulationVec[1]))	    
	    ## Product of probabilities:
	    binomPVec <- choose(nDiseased,d1Vec)*(pVec[1])^d1Vec * 
			(pVec[2])^(nDiseased-d1Vec) 
	    ## Compute significance of stratified sampling scheme:
	    out <- sum(binomPVec * sapply(d1Vec, function(x) computePValue(nPopulationVec[1],
			nSampleVec[1], x, sensitivity[1],specificity)) * 
			sapply((nDiseased-d1Vec), function(x) computePValue(nPopulationVec[2],
			nSampleVec[2], x, sensitivity[2],specificity)))		
	} else {
        ## Matrix of possible combinations of diseased elements in the groups:
	    dMx <- computeDiseaseMatrix(nDiseased, nPopulationVec)
		## Matrix with probabilities:
	    ## Powers of probabilities p_i^d_i:
	    pMx <- t(pVec^t(dMx))	
	    ## Hypergeometric probabilities:
	    pHypMx <- sapply(seq(along = dMx[1,]), function(jj){
    	    dMax <- min(c(nPopulationVec[jj], nDiseased))
	        dVec <- 0:dMax
	        pHypVec <- sapply(dVec, function(x) computePValue(nPopulationVec[jj],
				nSampleVec[jj], x, sensitivity[jj], specificity))
	        outVec <- pHypVec[(dMx[,jj]+1)]
        })
        ## Compute multinomial coefficients:
        multinomialCoeffVec <- sapply(seq(along = dMx[,1]), function(ii){
		    exp(lfactorial(sum(dMx[ii,])) - sum(lfactorial(dMx[ii,])))		
	    })
        ## Compute sensitivity of system:
        out <- sum(apply(pHypMx*pMx, prod, MARGIN = 1)*multinomialCoeffVec)
	}	
	## Return value:
    return(out)
}


