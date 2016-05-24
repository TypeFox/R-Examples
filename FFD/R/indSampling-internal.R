## Ian Kopacka
## 2010-07-16
##
## Function: indSampling.internal
## Package: FFD
## 
## For a fixed Population the number of herds to be tested, the
## expected number of animals to be tested and the costs are 
## computed for a sequence of herd sensitivities using "individual
## sampling".
##
## This internal function is called by 'indSamplingSummary.R'
##
## Input parameters:
##     survey.Data....Object of the class 'SurveyData', created by using
##                    the function 'surveyData.R'
##     herdSensVec....Numeric vector with values in (0,1). Herd sensitivities for
##                    which the number of herds to be tested, the expected total 
##                    number of animals to be tested and the expected costs are
##                    computed.
##    nSampleFixVec...Numeric vector containing some NAs (optional argument). 
##                    For risk groups for which the sample size is fixed 
##                    specify the sample size. For the risk groups for which
##                    the sample size should be computed set NA (order of the
##                    risk groups must be the same order as in 
##                    'survey.Data@riskValueData').
##    probVec.........Numeric vector. For those risk groups for which the 
##                    sample size should be computed sample probabilities must 
##                    be specified.
##                    The vector must have the same length as the number of 
##                    NA entries in nSampleFixVec or if nSampleFixVec is not 
##                    specified, probVec must have the same length as the 
##                    number of rows in survey.Data@riskValueData.
##
## Return value: Data frame with columns 'herdSensVec', 'nHerdsVec', 
##               'nAnimalsMeanVec' and 'expectedCostVec'.

## Internal auxilliary function:
indSampling.internal <- function(survey.Data, herdSensVec,
	nSampleFixVec = NULL, probVec = NULL){
    
    ## Check: Can risk based sampling be performed:
	if (!(is.null(nSampleFixVec) & is.null(probVec))){
	    if (dim(survey.Data@riskValueData)[1] == 0){
			warning (paste("'riskValueData' not specified in 'survey.Data'.\n", 
				"Performing sampling without risk groups"))
            nSampleFixVec <- NULL 
		    probVec <- NULL			
		}		
	}
	if (!(is.null(nSampleFixVec) & is.null(probVec))){
	    if (length(survey.Data@riskGroupVec) == 0){
			warning (paste("'riskGroupVec' not specified in 'survey.Data'.\n", 
				"Performing sampling without risk groups"))
            nSampleFixVec <- NULL 
		    probVec <- NULL			
		}		
	}

	## No risk groups:
	##################
	if (is.null(nSampleFixVec) & is.null(probVec)){	
		## Number of herds to be tested according to the
		## herd sensitivities:
		nHerdsVec <- sapply(herdSensVec, 
			function(x) computeOptimalSampleSize(nPopulation = length(survey.Data@nAnimalVec), 
				prevalence = survey.Data@designPrevalence, 
				alpha = survey.Data@alpha,
				sensitivity = x, specificity = 1, lookupTable = FALSE))
		
		## Only consider herd sensitivities for which the error level can be achieved:
		indVec <- is.finite(nHerdsVec)
		nHerdsVec <- nHerdsVec[indVec]
		herdSensVec <- herdSensVec[indVec]
		if (length(herdSensVec) == 0) return(NULL)
		
		## Number of Animals to be tested:
		meanSampleSize <- sapply(herdSensVec, 
				function(x){
					out <- computeSampleSizeInd(survey.Data = survey.Data, 
							herdSensitivity = x)
					return(out$nAnimalsMeanPerHerd)
				})        
		nAnimalsMeanVec <- meanSampleSize*nHerdsVec    
		
		## Number of herds to test per risk group (no risk groups):
		nHerdsMx <- matrix(numeric(), 0,0)
	} else {
	## With risk groups:
	####################		
		## Determine the number of herds in each risk group:
		riskValueDf <- survey.Data@riskValueData[,1:2]
		names(riskValueDf) <- c("riskGroup", "riskValues")
		riskValueDf$riskGroup <- as.character(riskValueDf$riskGroup)
		riskValueDf$id <- seq(along = riskValueDf[,1])
		riskGroupTab <- table(survey.Data@riskGroupVec)
		riskGroupDf <- data.frame(riskGroup = as.character(names(riskGroupTab)), 
				nPopulation = as.vector(riskGroupTab))
		riskValueDf <- merge(x = riskValueDf, y = riskGroupDf, by = "riskGroup",
				sort = FALSE)	
		riskValueDf <- riskValueDf[order(riskValueDf$id, decreasing = FALSE),]
		
		## Fixed sample sizes:
	    if (is.null(nSampleFixVec)){ 
	        nSampleFixVec <- rep(NA, length(riskValueDf$nPopulation))
		    class(nSampleFixVec) <- "numeric"
	    }
		
		## Compute the number of herds to be tested per risk group:
		nHerdsRG <- lapply(herdSensVec,
			function(x) computeOptimalSampleSizeRiskGroups(
				nPopulationVec = riskValueDf$nPopulation, 
				nRelRiskVec = riskValueDf$riskValues, 
				nSampleFixVec = nSampleFixVec,
				nSamplePropVec = probVec*riskValueDf$nPopulation[is.na(nSampleFixVec)],
				prevalence = survey.Data@designPrevalence, 
				alpha = survey.Data@alpha, 
				sensitivity = x, 
				specificity = 1))
		indVec <- Reduce(function(x,y) c(x,y), lapply(nHerdsRG, function(x) x[1]))
		indVec <- is.finite(indVec)
		nHerdsMx <- Reduce(function(x,y) rbind(x,y), nHerdsRG[indVec])
		row.names(nHerdsMx) <- NULL
		colnames(nHerdsMx) <- riskValueDf$riskGroup
		## Remove sample sizes for which confidence could not be reached:
		herdSensVec <- herdSensVec[indVec]
		if (length(herdSensVec) == 0) return(NULL)
		## Total number of herds to be tested:
		nHerdsVec <- rowSums(nHerdsMx)
		
		## Number of Animals to be tested:
		nAnimalsMeanVec <- sapply(seq(along = herdSensVec), 
			function(ii){
				tempList <- computeSampleSizeInd(survey.Data = survey.Data, 
			        herdSensitivity = herdSensVec[ii],
			        groupVec = survey.Data@riskGroupVec)
			    return(sum(tempList$nAnimalsMeanPerHerd[colnames(nHerdsMx)] * 
					nHerdsMx[ii,]))
		}) 		
	}
    ## Total Cost:
    if ((length(survey.Data@costHerd) > 0) & (length(survey.Data@costAnimal) > 0)){
        expectedCostVec <- nAnimalsMeanVec*survey.Data@costAnimal + 
            nHerdsVec*survey.Data@costHerd
    } else {
        expectedCostVec <- numeric(length(herdSensVec))*NA
    }
        
	## Return value:
	out <- list(herdSensVec = herdSensVec, 
        nHerdsVec = nHerdsVec,
		nHerdsPerRiskGroupMx = nHerdsMx,
		nAnimalsMeanVec = nAnimalsMeanVec, 
		expectedCostVec = expectedCostVec)
    return(out)    
}
