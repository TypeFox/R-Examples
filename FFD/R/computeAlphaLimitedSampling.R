## Ian Kopacka
## 2010-05-05
##
## Function: computeAlphaLimitedSampling
## 
## Compute the average alpha error for limited sampling.
##
## Input parameters:
##    stockSizeVector........Integer vector. Stock sizes of the herds.
##    sampleSizeLtd..........Integer. Sample size for limited sampling.
##    intraHerdPrevalence....Numeric between 0 and 1. Intra Herd prevalence.
##    diagSensitivity........Numeric between 0 and 1. Sensitivity of the 
##                           diagnostic test.
##    diagSpecificity........Numeric between 0 and 1. Specificity of the 
##                           diagnostic test.
##    groupVec...............Character vector. Optional parameter. If specified 
##                           it must have the same length as 'stockSizeVector'. 
##                           Defines the gouping of the data. Mean Alpha is 
##                           then returned for each group.
##
## Return value: List with slots
##    alphaDataFrame...Data frame. Variables "size" = vector of the unique 
##                     herd sizes in the population, "alpha" = alpha error
##                     made by limited sampling scheme for each herd size.
##    meanAlpha........Numeric between 0 and 1. Mean alpha error for the
##                     given population and sampling scheme. E.g. used
##                     for calculating the optimal sample size.
##
## Source: M. Ziller, T. Selhorst, J. Teuffert, M. Kramer and H. Schlueter, 
## "Analysis of sampling strategies to substantiate freedom from disease in 
## large areas", Prev. Vet. Med. 52 (2002), pp. 333--343.
##
## Calls:
##    computePValue.R
##
## Is called by:
##    -
##
computeAlphaLimitedSampling <- function(stockSizeVector, 
    sampleSizeLtd, intraHerdPrevalence, diagSensitivity,
    diagSpecificity = 1, groupVec = NULL){
    ## Vector of unique stock sizes:
    stockSizeUniqueVector <- unique(as.numeric(as.character(stockSizeVector)))
    ## Compute corresponding sample sizes (limited sampling):
    sampleSizesUniqueVector <- pmin.int(stockSizeUniqueVector, sampleSizeLtd)
    ## Compute the number of diseased individuals in the herds (minimum = 1):
    nDiseasedUniqueVector <- round(stockSizeUniqueVector * intraHerdPrevalence)
    nDiseasedUniqueVector <- pmax.int(nDiseasedUniqueVector, 1)
    ## Compute alpha errors for the different herd sizes:
    alphaUniqueVector <- sapply(seq(along = stockSizeUniqueVector), 
        function(ii){computePValue(nPopulation = 
        stockSizeUniqueVector[ii], nSample = sampleSizesUniqueVector[ii], 
        nDiseased = nDiseasedUniqueVector[ii], sensitivity = 
        diagSensitivity, specificity = diagSpecificity)})
    alphaSizesDataFrame <- data.frame(size = stockSizeUniqueVector, 
        alpha = alphaUniqueVector)
    alphaSizesDataFrame <- alphaSizesDataFrame[order(alphaSizesDataFrame$size,
		decreasing = FALSE),]    
    row.names(alphaSizesDataFrame) <- NULL

	## Compute the distribution of the herd sizes in the population 
    ## according to stockSizeVector, i.e., for each stock size 
    ## compute the probability of a herd being that particular size
    ## (= relative frequency).
    stockSizeTable <- table(stockSizeVector)
    stockSizeDistributionDataFrame <- data.frame(size = 
        as.numeric(as.character(names(stockSizeTable))), 
        distr = as.vector(stockSizeTable)/length(stockSizeVector))
    stockSizeDistributionDataFrame <- merge(x = stockSizeDistributionDataFrame, 
        y = alphaSizesDataFrame, by = "size", all = TRUE)
	meanAlpha <- sum(stockSizeDistributionDataFrame$distr * 
        stockSizeDistributionDataFrame$alpha)
	meanAlphaRiskGroups <- NULL
	
    ## Grouping variable specified:
	if (!is.null(groupVec)){		
	    splitList <- split(x = stockSizeVector, f = groupVec)
		meanAlphaList <- lapply(splitList, function(nAnimalVec){
		    stockSizeTable <- table(nAnimalVec)
			stockSizeDistributionDataFrame <- data.frame(size = 
                as.numeric(as.character(names(stockSizeTable))), 
                distr = as.vector(stockSizeTable)/length(nAnimalVec))
            stockSizeDistributionDataFrame <- merge(x = stockSizeDistributionDataFrame, 
                y = alphaSizesDataFrame, by = "size", all.x = TRUE, all.y = FALSE)
	        meanGroupAlpha <- sum(stockSizeDistributionDataFrame$distr * 
                stockSizeDistributionDataFrame$alpha)
		})
        meanAlphaRiskGroups <- Reduce(function(x,y) c(x,y), meanAlphaList)
	    names(meanAlphaRiskGroups) <- names(meanAlphaList)	
	}
    ## Berechne gemittelten alpha-Fehler:
    out <- list(alphaDataFrame = alphaSizesDataFrame, meanAlpha = meanAlpha,
		meanAlphaRiskGroups = meanAlphaRiskGroups)    
    return(out)    
}
