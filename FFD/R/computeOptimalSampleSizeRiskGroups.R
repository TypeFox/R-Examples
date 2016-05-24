## Ian Kopacka
## 2011-07-05
##
## Function: computeOptimalSampleSizeRiskGroups
## 
## Computes the optimal sample size (for each risk group) for a survey to 
## substantiate freedom from disease for a population stratified into risk groups. 
## The optimal sample size is the smallest sample size that produces 
## an alpha-error less than or equal to a prediscribed value for alpha. 
## The population is considered as diseased if at least one individual has a 
## positive test result. The sample size is computed using a bisection method.
##
## The sample size can be fixed for a subset of the risk groups via the input 
## parameter 'nSampleFixVec' (vector containing sample sizes for the risk groups
## with fixed values and NA for the risk groups for which the sample size is to 
## be computed). For those risk groups for which the sample size is to be
## computed a vector specifying the proportional distribution among the risk 
## groups ('nSamplePropVec') needs to be specified.
##
## Example: We have 3 risk groups. For the 2nd risk group we want 20 farms
## to be sampled. For the other  risk groups we specify that the sample size 
## for risk group 1 should be double the sample size of risk group 3. We
## then set:
##     nSampleFixVec <- c(NA, 20, NA)
##     nSamplePropVec <- c(2,1)
##
## Input parameters:
##     nPopulationVec...Integer vector. Population sizes of the risk groups.
##     nRelRiskVec......Numeric vector. (Relative) infection risks of the 
##                      risk groups.
##     nSampleFixVec.......Numeric vector containing NAs (optional argument). 
##                      For risk groups for which the sample size is fixed 
##                      specify the sample size. For the risk groups for which
##                      the sample size should be computed set NA (order of the
##                      risk groups must be the same order as in 'nPopulationVec'
##                      and 'nRelRiskVec') .
##     nSamplePropVec...Numeric vector. For those risk groups for which the 
##                      sample size should be computed a proportional 
##                      distribution of the overall sample size must be specified.
##                      The vector must have the same length as the number of 
##                      NA entries in nSampleFixVec or if nSampleFixVec is not 
##                      specified, nSamplePropVec must have the same length as
##                      nPopulationVec. 
##     prevalence.......Numeric between 0 and 1. Prevalence of the disease 
##                      under the null hypothesis.
##     alpha............Numeric between 0 and 1. Type one error for the 
##                      statistical test.
##     sensitivity......Numeric between 0 and 1. Sensitivity of test (diagnostic 
##                      test for one stage sampling, herd test for two stage 
##                      sampling).
##     specificity......Numeric between 0 and 1. Specificity of test (diagnostic 
##                      test for one stage sampling, herd test for two stage 
##                      sampling).
##     lookupTable......Logical. TRUE if a lookup table of sample sizes for
##                      individual sampling (see, Ziller et al., 2002) should
##                      be produced. FALSE if the sample size is desired
##                      for a fixed population size (default).
##
## Calls: 
##    computePValueRiskGroups.R
##
## Is called by:
##    -
##
computeOptimalSampleSizeRiskGroups <- function(nPopulationVec, nRelRiskVec, 
	nSampleFixVec = NULL, nSamplePropVec = NULL, prevalence, alpha = 0.05, 
    sensitivity = 1, specificity = 1){
    ## Fixed sample sizes:
	if (is.null(nSampleFixVec)) nSampleFixVec <- rep(NA, length(nPopulationVec))
	
	## Error check:
	if (length(nRelRiskVec) != length(nPopulationVec)){
		stop ("'nPopulationVec' and 'nRelRiskVec' must have the same length.")		
	}
	if (length(nSampleFixVec) != length(nPopulationVec)){
		stop ("'nSampleFixVec' and 'nRelRiskVec' must have the same length.")	
	}
	if ((sum(is.na(nSampleFixVec)) > 1) & (length(nSamplePropVec) != sum(is.na(nSampleFixVec)))){
		stop ("'nSamplePropVec' does not have the correct length.")		
	}	
	if ((sum(is.na(nSampleFixVec)) == 1) & (!is.null(nSamplePropVec))){
		warning ("Argument 'nSamplePropVec' is ignored.")		
	}
	if (min(nSamplePropVec) <= 0){
		stop ("'nSamplePropVec' must have positive values.")		
	}
	if (sum(is.na(nSampleFixVec)) == 0){
		stop ("'nSampleFixVec' must contain NA values.")		
	}

    ## Calculate the number of diseased individuals in the population:
    nDiseased <- max(round(sum(nPopulationVec)*prevalence),1)
	## Proportion for sample sizes:
	if (sum(is.na(nSampleFixVec)) > 1){
	    nSamplePropVec <- nSamplePropVec/sum(nSamplePropVec)
	} else {
		nSamplePropVec <- 1
	}	
    
	# Initialize the parameters for the bisection method:
    if (sum(nSampleFixVec, na.rm = TRUE) > 0){
	    lowerBoundSampleSize <- 0
	} else {
		lowerBoundSampleSize <- 1
	}	
    upperBoundSampleSize <- min(floor(nPopulationVec[is.na(nSampleFixVec)] / 
		nSamplePropVec))
    nSampleArgVec <- nSampleFixVec
    ## Compute alpha errors:
	probabilityLowerUpper <- sapply(c(lowerBoundSampleSize, 
        upperBoundSampleSize), function(sampleSize){ 	    
		nSampleArgVec[is.na(nSampleFixVec)] <- 
			roundConstantSum(sampleSize*nSamplePropVec, output = 0)
	    out <- computePValueRiskGroups(nPopulationVec = nPopulationVec, 
			nSampleVec = nSampleArgVec, nRelRiskVec = nRelRiskVec,
	        nDiseased = nDiseased, sensitivity = sensitivity, 
			specificity = specificity)
	    return(out)	
    })
    ## Check if upper bound satisfies the desired accuracy:
    if (probabilityLowerUpper[2] > alpha) return(Inf)
    ## Check if lower bound satisfies the desired accuracy:
    if (probabilityLowerUpper[1] <= alpha){
		nSampleOutVec <- nSampleFixVec
	    nSampleOutVec[is.na(nSampleFixVec)] <- 
		    roundConstantSum(lowerBoundSampleSize*nSamplePropVec, output = 0)
	    return(nSampleOutVec)
	} 
	## Bisection method: Due to the discrete nature of the variable
    ## (= optimal sample size) the bisection method is executed
    ## until the width of the search interval falls below a certain
    ## threshold. The probabilities for the remaining cases are 
    ## explicitly computed and the optimal sample size is chosen:
    minIntervalWidth <- 4 #10
    while (upperBoundSampleSize - lowerBoundSampleSize > minIntervalWidth){
        midSampleSize <- round((lowerBoundSampleSize + upperBoundSampleSize)/2)
        nSampleArgVec[is.na(nSampleFixVec)] <- 
			roundConstantSum(midSampleSize*nSamplePropVec, output = 0)
	    probabilityMid <- computePValueRiskGroups(nPopulationVec = nPopulationVec, 
			nSampleVec = nSampleArgVec, nRelRiskVec = nRelRiskVec,
	        nDiseased = nDiseased, sensitivity = sensitivity, 
			specificity = specificity)
        if (probabilityMid > alpha){
            lowerBoundSampleSize <- midSampleSize
        } else {
            upperBoundSampleSize <- midSampleSize
        }    
    }
	## Explicitly compute the probabilities for the remaining cases
    ## and compare the cases. Return the smallest sample size that
    ## has a probability of returning no test positives smaller than
    ## or equal to alpha:
    sampleSizeVector <- lowerBoundSampleSize : upperBoundSampleSize
    probabilityVector <- sapply(sampleSizeVector, 
        function(sampleSize){ 
		nSampleArgVec <- nSampleFixVec
		nSampleArgVec[is.na(nSampleFixVec)] <- 
			roundConstantSum(sampleSize*nSamplePropVec, output = 0)
	    out <- computePValueRiskGroups(nPopulationVec = nPopulationVec, 
			nSampleVec = nSampleArgVec, nRelRiskVec = nRelRiskVec,
	        nDiseased = nDiseased, sensitivity = sensitivity, 
			specificity = specificity)
	    return(out)	
    })
    indexProbablityAlpha <- which(probabilityVector <= alpha)
    if (length(indexProbablityAlpha) > 0){
        sampleSizeOpt <- sampleSizeVector[min(indexProbablityAlpha)]
    } else {
        sampleSizeOpt <- upperBoundSampleSize
    }
	## Return value = optimal sample size:
	nSampleOutVec <- nSampleFixVec
	nSampleOutVec[is.na(nSampleFixVec)] <- 
		roundConstantSum(sampleSizeOpt*nSamplePropVec, output = 0)
	return(nSampleOutVec)
}
