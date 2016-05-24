## Ian Kopacka
## 2011-07-07
##
## Function: indSampling
## 
## Constructor for objects of the class 'IndSampling'.
##
## Package: FFD
##
## Input parameters:
##    survey.Data.......Object of the class 'SurveyData', created by using
##                      the function 'surveyData.R'
##    herdSensitivity...Numeric between 0 and 1. Desired herd sensitivity.
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
## Return value: object of the class 'IndSampling'.

indSampling <- function(survey.Data, herdSensitivity, nSampleFixVec = NULL, 
	probVec = NULL){
    # Error check:
    if ((herdSensitivity >= 1) | (herdSensitivity <= 0)){ 
		stop("[indSampling]: herdSensitivity must be a value in (0,1).\n")
	}
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
	if (is.null(nSampleFixVec) & is.null(probVec)){		
		## Number of herds to test:
		nHerds <- computeOptimalSampleSize(nPopulation = length(survey.Data@nAnimalVec), 
				prevalence = survey.Data@designPrevalence, alpha = survey.Data@alpha,
				sensitivity = herdSensitivity, specificity = 1, lookupTable = FALSE)
		
		## Lookup table for the number of animals to test per herd depending on the 
		## herd size and expected total number of animals to test:
		tempList <- computeSampleSizeInd(survey.Data = survey.Data, 
				herdSensitivity = herdSensitivity)
		nAnimalsMean <- tempList$nAnimalsMeanPerHerd*nHerds
		lookupTable <- tempList$lookupTable
		## Number of herds per risk group (no risk groups used):
		nHerdsPerRiskGroup <- numeric()
		nSampleFixVec <- numeric()
		probVec <- numeric()
	} else {
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
	    ## Compute the number of herds to be tested:
		nHerdsRG <- computeOptimalSampleSizeRiskGroups(
			nPopulationVec = riskValueDf$nPopulation, 
			nRelRiskVec = riskValueDf$riskValues, 
			nSampleFixVec = nSampleFixVec,
			nSamplePropVec = probVec*riskValueDf$nPopulation[is.na(nSampleFixVec)],
			prevalence = survey.Data@designPrevalence, 
			alpha = survey.Data@alpha, 
			sensitivity = herdSensitivity, 
			specificity = 1)
		if (nHerdsRG[1] == Inf){ 
			stop(paste("Desired alpha cannot be reached", 
				"with specified 'nSampleFixVec' and 'probVec'."))
		}
		#survey.Data@riskValueData$nHerds <- nHerdsRG
		nHerds <- sum(nHerdsRG)
		## Lookup table for the number of animals to test per herd depending on the 
		## herd size and expected total number of animals to test:
		tempList <- computeSampleSizeInd(survey.Data = survey.Data, 
			herdSensitivity = herdSensitivity,
			groupVec = survey.Data@riskGroupVec)
	    nAnimalsMeanDf <- data.frame(riskGroup = names(tempList$nAnimalsMeanPerHerd),
			nAnimalsMean = as.vector(tempList$nAnimalsMeanPerHerd))
	    nAnimalsMeanDf <- merge(x = riskValueDf, y = nAnimalsMeanDf, 
			by = "riskGroup", sort = FALSE)
	    nAnimalsMeanDf <- nAnimalsMeanDf[order(nAnimalsMeanDf$id, decreasing = FALSE),]
		nAnimalsMean <- sum(nHerdsRG*nAnimalsMeanDf$nAnimalsMean)
		lookupTable <- tempList$lookupTable
		## Number of herds per risk group (no risk groups used):
		nHerdsPerRiskGroup <- nHerdsRG
		names(nHerdsPerRiskGroup) <- riskValueDf$riskGroup			
	}
	
    ## Expected cost:
    expectedCost <- nHerds*survey.Data@costHerd + nAnimalsMean*survey.Data@costAnimal
    
    ## Create object of class 'IndSampling':
    out <- new("IndSampling", surveyData = survey.Data,
        herdSensitivity = herdSensitivity,
        nHerds = nHerds, 
		nHerdsPerRiskGroup = nHerdsPerRiskGroup,
		nSampleFixVec = nSampleFixVec,
		probVec = probVec,
        nAnimalsMean = nAnimalsMean,
        expectedCost = expectedCost,
        lookupTable = lookupTable)
    return(out)    
}
