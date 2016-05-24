## Ian Kopacka
## 2011-07-05
##
## Function: ltdSampling
## 
## Constructor for objects of the class 'LtdSampling'. If values for 
## 'nSampleFixVec' and/or 'probVec' are specified, sampling with stratification
## by risk groups is performed, otherwise herds are sampled without the 
## use of risk groups.
##
## Package: FFD
##
## Input parameters:
##    survey.Data.......Object of the class 'SurveyData', created by using
##                      the function 'surveyData.R'
##    sampleSizeLtd.....Positive integer. Pre-fixed number of animals to be 
##                      tested per holding, irrespective of the herd size 
##                      (if the herd contains fewer animals then the entire
##                      herd needs to be tested).
##    nSampleFixVec.....Numeric vector containing some NAs (optional argument). 
##                      For risk groups for which the sample size is fixed 
##                      specify the sample size. For the risk groups for which
##                      the sample size should be computed set NA (order of the
##                      risk groups must be the same order as in 
##                      'survey.Data@riskValueData').
##    probVec...........Numeric vector. For those risk groups for which the 
##                      sample size should be computed sample probabilities must 
##                      be specified.
##                      The vector must have the same length as the number of 
##                      NA entries in nSampleFixVec or if nSampleFixVec is not 
##                      specified, probVec must have the same length as the 
##                      number of rows in survey.Data@riskValueData.
##
## Return value: object of the class 'LtdSampling'.

ltdSampling <- function(survey.Data, sampleSizeLtd, nSampleFixVec = NULL, 
		probVec = NULL){
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
		## Compute the mean herd sensitivity:
	    alphaList <- computeAlphaLimitedSampling(stockSizeVector = survey.Data@nAnimalVec, 
			sampleSizeLtd = sampleSizeLtd, 
			intraHerdPrevalence = survey.Data@intraHerdPrevalence, 
			diagSensitivity = survey.Data@diagSensitivity, 
			diagSpecificity = 1)		
		## Compute the number of herds to be tested:
		nHerds <- computeOptimalSampleSize(nPopulation = length(survey.Data@nAnimalVec), 
				prevalence = survey.Data@designPrevalence, 
				alpha = survey.Data@alpha, 
				sensitivity = (1 - alphaList$meanAlpha), 
				specificity = 1, lookupTable = FALSE)
		if (nHerds == Inf){ 
			stop(paste("Desired alpha cannot be reached", 
							"even when entire population is tested."))
		}
		## Compute mean overall number of Animals to be tested:
		nAnimalsMean <- nHerds*mean(pmin(survey.Data@nAnimalVec,sampleSizeLtd))
		## Number of herds per risk group (no risk groups used):
		nHerdsPerRiskGroup <- numeric()
		meanHerdSensPerRG <- numeric()
		nSampleFixVec <- numeric()
		probVec <- numeric()
		
	} else {	
	## With risk groups:
	####################
	    ## Compute the mean herd sensitivity:
	    alphaList <- computeAlphaLimitedSampling(stockSizeVector = survey.Data@nAnimalVec, 
			sampleSizeLtd = sampleSizeLtd, 
			intraHerdPrevalence = survey.Data@intraHerdPrevalence, 
			diagSensitivity = survey.Data@diagSensitivity, 
			diagSpecificity = 1,
			groupVec = survey.Data@riskGroupVec)
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
		riskValueDf$meanHerdSens <- as.vector(1 - 
			alphaList$meanAlphaRiskGroups[riskValueDf$riskGroup])
		## Fixed sample sizes:
	    if (is.null(nSampleFixVec)){ 
			nSampleFixVec <- rep(NA, length(riskValueDf$nPopulation))
			class(nSampleFixVec) <- "numeric"
		}
	    ## Compute the number of herds to be tested:
		nHerdsRG <- computeOptimalSampleSizeRiskGroups(
			nPopulationVec = riskValueDf$nPopulation, 
			nRelRiskVec = riskValueDf$riskValues, 
			nSampleFixVec = nSampleFixVec,
			nSamplePropVec = probVec*riskValueDf$nPopulation[is.na(nSampleFixVec)],
			prevalence = survey.Data@designPrevalence, 
			alpha = survey.Data@alpha, 
			#sensitivity = (1 - alphaList$meanAlpha),
			sensitivity = riskValueDf$meanHerdSens, 
			specificity = 1)
		if (nHerdsRG[1] == Inf){ 
			stop(paste("Desired alpha cannot be reached", 
				"with specified 'nSampleFixVec' and 'probVec'."))
		}
		#survey.Data@riskValueData$nHerds <- nHerdsRG
		nHerds <- sum(nHerdsRG)
		## Compute mean overall number of Animals to be tested:
		nAnimalTab <- tapply(X = survey.Data@nAnimalVec, INDEX = survey.Data@riskGroupVec,
			FUN = function(x) mean(pmin(x,sampleSizeLtd)))
	    nAnimalDf <- data.frame(riskGroup = as.character(names(nAnimalTab)), 
			nAnimalsMean = as.vector(nAnimalTab))
		riskValueDf <- merge(x = riskValueDf, y = nAnimalDf, by = "riskGroup",
		    sort = FALSE)	
	    riskValueDf <- riskValueDf[order(riskValueDf$id, decreasing = FALSE),]
		nAnimalsMean <- sum(nHerdsRG*riskValueDf$nAnimalsMean)	
		## Number of herds per risk group (no risk groups used):
		nHerdsPerRiskGroup <- nHerdsRG
		names(nHerdsPerRiskGroup) <- riskValueDf$riskGroup		
		meanHerdSensPerRG <- riskValueDf$meanHerdSens
		names(meanHerdSensPerRG) <- riskValueDf$riskGroup		
		
	}
	## Expected cost of the survey:
	expectedCost <- nHerds*survey.Data@costHerd + 
		nAnimalsMean*survey.Data@costAnimal
	## Create object of class 'LtdSampling':
	out <- new("LtdSampling", surveyData = survey.Data,
			sampleSizeLtd = sampleSizeLtd,
			meanHerdSensitivity = (1 - alphaList$meanAlpha),
			meanHerdSensPerRG = meanHerdSensPerRG,
			nHerds = nHerds, 
			nHerdsPerRiskGroup = nHerdsPerRiskGroup,
			nSampleFixVec = nSampleFixVec,
			probVec = probVec,
			nAnimalsMean = nAnimalsMean,
			expectedCost = expectedCost)
	return(out)     
}
