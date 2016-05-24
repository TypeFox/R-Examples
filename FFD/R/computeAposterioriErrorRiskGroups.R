## Ian Kopacka
## 2011-07-08
##
## Function: computeAposterioriErrorRiskGroups
## 
## Computes the a posteriori error, i.e., the probability of finding no 
## testpositives in the sample, given the disease is present at the design 
## prevalence. The probability is computed for a population stratified
## by a risk factor with two or more possible values, for which the relative
## risks of being infected are known.
##
## Input parameters:
##   alphaErrorVector...Numeric vector. Alpha-error (= 1-herd sensitivity)
##                      of the herds in the sample.
##   groupVec...........Character vector. Risk group to which each of the 
##                      herds in the sample belongs. Must have the same length 
##                      (and order) as 'alphaErrorVector'.
##   nPopulationVec.....Integer vector. Population sizes of the risk groups.
##   nRelRiskVec........Numeric vector. (Relative) infection risks of the 
##                      risk groups.
##   nDiseased..........Integer. Number of diseased in the total population 
##                      (according to the null hypothesis).
##   method............."exact" for exact error, "approx" for approximation
##                      recommended for nDiseased > 7
##
## Sources: 
## A.R. Cameron, F.C. Baldock, "A new probability fomula for 
## surveys to substantiate freedom from disease", Prev. Vet. Med. 34
## (1998), pp. 1 - 17.
## P.A.J. Martin, A.R. Cameron, M. Greiner, "Demonstrating freedom from disease
## using multiple complex data sources. !: A new methodology based on scenario
## trees", Prev. Vet. Med. 79 (2007), pp. 71 - 97.
##
## Calls:
##    
##
## Is called by:
##

computeAposterioriErrorRiskGroups <- function(alphaErrorVector, groupVec,
	groupLevels, nPopulationVec, nRelRiskVec, nDiseased, method = "default"){
    
    ## Determine sample method:
    if (!(method %in% c("exact", "approx", "approximation", "default"))){
        stop("Undefined value for argument 'method'; must be either 'exact' or 'approx'.")
    } 
    if (method == "default"){
        if (nDiseased > 3){
            method <- "approx"
        } else {
            method <- "exact"
        }
    }
	
	## Error check:
    nSample <- length(alphaErrorVector)  # sample size
	nPopulation <- sum(nPopulationVec)
    if(nPopulation < nSample) stop("Sample size is larger than population size.")
    if(nPopulation < nDiseased) stop("Number of diseased is larger than population size.")
    if(length(alphaErrorVector) != length(groupVec)) stop("'alphaErrorVector' must have the same length as 'groupVec'")
    if(length(nPopulationVec) != length(nRelRiskVec)) stop("'nPopulationVec' must have the same length as 'nRelRiskVec'")
    if((min(nPopulationVec) <= 0) | any (nPopulationVec - as.integer(nPopulationVec) != 0)) stop("'nPopulationVec' must be a positive integer.")
	if((min(nDiseased) <= 0) | any (nDiseased - as.integer(nDiseased) != 0)) stop("'nDiseased' must be a positive integer.")
	
    ## Factor for normalization:
	theta <- sum(nPopulationVec*nRelRiskVec)
	## Vector of probabilities p_i:
	pVec <- nRelRiskVec*nPopulationVec/theta 

	## Two Risk groups:
	if (length(nPopulationVec) == 2){
		## Possible number of diseased in RG 1:
		d1Vec <- max(c(0,nDiseased-nPopulationVec[2])):
			min(c(nDiseased,nPopulationVec[1]))	    
	    ## Product of probabilities:
	    binomPVec <- choose(nDiseased,d1Vec)*(pVec[1])^d1Vec * 
			(pVec[2])^(nDiseased-d1Vec) 
	    ## Compute a posteriori error:
	    out <- sum(binomPVec * sapply(d1Vec, function(x){
		    alpha1 <- computeAposterioriError(alphaErrorVector = 
				alphaErrorVector[groupVec == groupLevels[1]], 
		        nPopulation = nPopulationVec[1], 
				nDiseased = x, method = method)	
		    alpha2 <- computeAposterioriError(alphaErrorVector = 
				alphaErrorVector[groupVec == groupLevels[2]], 
		        nPopulation = nPopulationVec[2], 
				nDiseased = nDiseased - x, method = method)
		    return(alpha1*alpha2)
			} 
        ))		
	} else {
#        ## Matrix of possible combinations of diseased elements in the groups:
#	    dMx <- computeDiseaseMatrix(nDiseased, nPopulationVec)
#		## Matrix with probabilities:
#	    ## Powers of probabilities p_i^d_i:
#	    pMx <- t(pVec^t(dMx))	
#	    ## Hypergeometric probabilities:
#	    pHypMx <- sapply(seq(along = dMx[1,]), function(jj){
#    	    dMax <- min(c(nPopulationVec[jj], nDiseased))
#	        dVec <- 0:dMax
#	        pHypVec <- sapply(dVec, function(x) computePValue(nPopulationVec[jj],
#				nSampleVec[jj], x, sensitivity[jj], specificity))
#	        outVec <- pHypVec[(dMx[,jj]+1)]
#        })
#        ## Compute multinomial coefficients:
#        multinomialCoeffVec <- sapply(seq(along = dMx[,1]), function(ii){
#		    exp(lfactorial(sum(dMx[ii,])) - sum(lfactorial(dMx[ii,])))		
#	    })
#        ## Compute sensitivity of system:
#        out <- sum(apply(pHypMx*pMx, prod, MARGIN = 1)*multinomialCoeffVec)
	}	
	## Return value:
    return(out)
	
	
}


