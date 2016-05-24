## Ian Kopacka
## 2010-07-12
##
## Function: ltdSampling.internal
## Package: FFD
## 
## For a fixed Population the number of herds to be tested, the
## expected number of animals to be tested and the costs are 
## computed for a sequence of sampling sizes using "limited
## sampling".
##
## This internal function is called by 'ltdSamplingSummary.R'
##
## Input parameters:
##     survey.Data........Object of the class 'SurveyData', created by using
##                        the function 'surveyData.R'
##     sampleSizeLtdVec...Numeric (integer) vector. Sample limit vector for
##                        which the mean herd sensitivity, the number of
##                        herds to be tested, the expected total number of
##                        animals to be tested and the expected costs are
##                        computed.
##    nSampleFixVec.......Numeric vector containing some NAs (optional argument). 
##                        For risk groups for which the sample size is fixed 
##                        specify the sample size. For the risk groups for which
##                        the sample size should be computed set NA (order of the
##                        risk groups must be the same order as in 
##                        'survey.Data@riskValueData').
##    probVec.............Numeric vector. For those risk groups for which the 
##                        sample size should be computed sample probabilities must 
##                        be specified.
##                        The vector must have the same length as the number of 
##                        NA entries in nSampleFixVec or if nSampleFixVec is not 
##                        specified, probVec must have the same length as the 
##                        number of rows in survey.Data@riskValueData.
##
## Return value: Data frame with columns 'sampleSizeLtdVec', 'meanHerdSensVec', 
##               'nHerdsVec', 'nAnimalsMeanVec' and 'expectedCostVec'.

## Internal auxilliary function:
ltdSampling.internal <- function(survey.Data, sampleSizeLtdVec, 
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
		## Compute the herd sensitivities:
	    meanAlphaVec <- sapply(sampleSizeLtdVec, function(x){
		    out <- computeAlphaLimitedSampling(stockSizeVector = survey.Data@nAnimalVec,
				sampleSizeLtd = x, intraHerdPrevalence = survey.Data@intraHerdPrevalence,
				diagSensitivity = survey.Data@diagSensitivity,
				diagSpecificity = 1)
			return(out$meanAlpha)
		})
	    meanHerdSensVec <- 1 - meanAlphaVec  
	
		## Number of herds to be tested according to the
		## herd sensitivities:
		nHerdsVec <- sapply(meanHerdSensVec, 
				function(x) computeOptimalSampleSize(nPopulation = length(survey.Data@nAnimalVec), 
							prevalence = survey.Data@designPrevalence, alpha = survey.Data@alpha,
							sensitivity = x, specificity = 1, lookupTable = FALSE))
		## Remove sample sizes for which confidence could not be reached:
		indVec <- is.finite(nHerdsVec)
		nHerdsVec <- nHerdsVec[indVec]
		meanHerdSensVec <- meanHerdSensVec[indVec]
		sampleSizeLtdVec <- sampleSizeLtdVec[indVec]
		if (length(sampleSizeLtdVec) == 0) return(NULL)
		
		## Mean total number of animals to be tested:
		nAnimalsMeanVec <- sapply(sampleSizeLtdVec, function(x) mean(pmin(x,survey.Data@nAnimalVec)))*nHerdsVec
		
		## Number of herds to test per risk group (no risk groups):
		nHerdsMx <- matrix(numeric(), 0,0)
		## Mean herd sensitivity per risk group (no risk groups):
		meanHerdSensPerRGMx <- matrix(numeric(), 0,0)
		
	} else {
	## With risk groups:
	####################	
	    meanAlphaMx <- t(sapply(sampleSizeLtdVec, function(x){
		    outList <- computeAlphaLimitedSampling(stockSizeVector = survey.Data@nAnimalVec,
				sampleSizeLtd = x, intraHerdPrevalence = survey.Data@intraHerdPrevalence,
				diagSensitivity = survey.Data@diagSensitivity,
				diagSpecificity = 1, groupVec = survey.Data@riskGroupVec)
		    outVec <- c(outList$meanAlpha, outList$meanAlphaRiskGroups)
			names(outVec) <- c("meanAlpha", names(outList$meanAlphaRiskGroups))
			return(outVec)
		}))
        meanHerdSensMx <- 1 - meanAlphaMx
	    meanHerdSensVec <- meanHerdSensMx[,"meanAlpha"]  
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
		nHerdsRG <- lapply(seq(along = meanHerdSensVec),
				function(ii) computeOptimalSampleSizeRiskGroups(
							nPopulationVec = riskValueDf$nPopulation, 
							nRelRiskVec = riskValueDf$riskValues, 
							nSampleFixVec = nSampleFixVec,
							nSamplePropVec = probVec*riskValueDf$nPopulation[is.na(nSampleFixVec)],
							prevalence = survey.Data@designPrevalence, 
							alpha = survey.Data@alpha, 
							sensitivity = as.vector(meanHerdSensMx[ii,riskValueDf$riskGroup]), 
							specificity = 1))
		indVec <- Reduce(function(x,y) c(x,y), lapply(nHerdsRG, function(x) x[1]))
		indVec <- is.finite(indVec)
		nHerdsMx <- Reduce(function(x,y) rbind(x,y), nHerdsRG[indVec])
		row.names(nHerdsMx) <- NULL
		colnames(nHerdsMx) <- riskValueDf$riskGroup
		## Remove sample sizes for which confidence could not be reached:
		sampleSizeLtdVec <- sampleSizeLtdVec[indVec]
		meanHerdSensVec <- meanHerdSensVec[indVec]
		meanHerdSensMx <- meanHerdSensMx[indVec,]
		meanHerdSensPerRGMx <- meanHerdSensMx[,riskValueDf$riskGroup]
		if (length(sampleSizeLtdVec) == 0) return(NULL)
		## Total number of herds to be tested:
		nHerdsVec <- rowSums(nHerdsMx)
		
		## Mean total number of animals to be tested:
		nAnimalsMeanVec <- sapply(seq(along = sampleSizeLtdVec),  
			function(ii){ 
				nAnimalTab <- tapply(X = survey.Data@nAnimalVec, 
					INDEX = survey.Data@riskGroupVec,
				    FUN = function(x) mean(pmin(x,sampleSizeLtdVec[ii])))
			    nAnimalDf <- data.frame(riskGroup = as.character(names(nAnimalTab)), 
				    nAnimalsMean = as.vector(nAnimalTab))
			    tempDf <- merge(x = riskValueDf, y = nAnimalDf, by = "riskGroup",
				    sort = FALSE)
			    tempDf <- tempDf[order(tempDf$id, decreasing = FALSE),]
				nAnimalsOut <- sum(nHerdsMx[ii,]*tempDf$nAnimalsMean)				
#			    nAnimalsOut <- tempDf$nAnimalsMean
			}	
		)					
	}
	## Total Cost:
	if ((length(survey.Data@costHerd) > 0) & (length(survey.Data@costAnimal) > 0)){
		expectedCostVec <- nAnimalsMeanVec*survey.Data@costAnimal + 
				nHerdsVec*survey.Data@costHerd
	} else {
		expectedCostVec <- numeric(length(sampleSizeLtdVec))*NA
	}
	
	## Return value:
	out <- list(sampleSizeLtdVec = sampleSizeLtdVec, 
			meanHerdSensVec = meanHerdSensVec, 
			meanHerdSensPerRGMx = meanHerdSensPerRGMx,
			nHerdsVec = nHerdsVec,
			nHerdsPerRiskGroupMx = nHerdsMx,
			nAnimalsMeanVec = nAnimalsMeanVec, 
			expectedCostVec = expectedCostVec)
	return(out)    
}
