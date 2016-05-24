## Ian Kopacka
## 2010-07-16
##
## Function: computeSampleSizeInd
## Package: FFD
## 
## For a fixed survey parameters and individual sampling the 
## function computes the lookup table for the number of animals 
## to test per herd depending on the herd size, as well as the
## average number of animals to test per herd according to the
## herd size distribution in the population.
##
## This internal function is called by 'indSampling-internal()'
## and 'indSampling()'
##
## Input parameters:
##    survey.Data.......Object of the class 'SurveyData', created by using
##                      the function 'surveyData.R'
##    herdSensitivity...Numeric between 0 and 1. Desired herd sensitivity.
##
## Return value: List with two list elements:
##
##    lookupTable...........Matrix with columns 'N_lower', 'N_upper', 
##                          'sampleSize'. Lookup table for the number of animals
##                          to test, depending on the herd size.
##
##    nAnimalsMeanPerHerd...Numeric. Mean number of animals to test
##                          per herd.

computeSampleSizeInd <- function(survey.Data, herdSensitivity, groupVec = NULL){
    ## Number of animals:
    #####################
    ## Compute lookup table:
    lookupTable <- computeOptimalSampleSize(nPopulation = max(survey.Data@nAnimalVec), 
        prevalence = survey.Data@intraHerdPrevalence, alpha = (1-herdSensitivity),
        sensitivity = survey.Data@diagSensitivity, specificity = 1, lookupTable = TRUE)
    nAnimalLookup <- as.data.frame(lookupTable)
    
    ## Merge with animal data:
    nAnimalLookup$interval <- paste("(", nAnimalLookup$N_lower-1, ",",
        nAnimalLookup$N_upper, "]", sep = "")
    breaks <- c(nAnimalLookup$N_lower[1]-1, nAnimalLookup$N_upper)     
	
	## No grouping Variable:
	if (is.null(groupVec)){
		nAnimalTable <- table(survey.Data@nAnimalVec)
		nAnimalDataFrame <- data.frame(nAnimal = as.numeric(as.character(names(nAnimalTable))),
				freq = as.vector(nAnimalTable), interval = cut(x = as.numeric(as.character(names(nAnimalTable))), 
						breaks = breaks, dig.lab = 10))        
		nAnimalDataFrame <- merge(x = nAnimalDataFrame, 
				y = subset(nAnimalLookup, select = c("interval", "sampleSize")),
				by = "interval", all.x = TRUE, all.y = FALSE)    
		## Mean number of animals to be tested per holding:
		nAnimalsMeanPerHerd <- sum(with(nAnimalDataFrame, freq*sampleSize))/length(survey.Data@nAnimalVec)
	} else {
	## Grouping variable specified:
	    splitList <- split(x = survey.Data@nAnimalVec, f = groupVec)
    	nAnimalsMeanPerHerdList <- lapply(splitList, function(nAnimalVec){
	        nAnimalTable <- table(nAnimalVec)			
		    nAnimalDataFrame <- data.frame(nAnimal = as.numeric(as.character(names(nAnimalTable))),
	    		freq = as.vector(nAnimalTable), 
			    interval = cut(x = as.numeric(as.character(names(nAnimalTable))), 
			    breaks = breaks, dig.lab = 10)) 
            nAnimalDataFrame <- merge(x = nAnimalDataFrame, 
			    y = subset(nAnimalLookup, select = c("interval", "sampleSize")),
			    by = "interval", all.x = TRUE, all.y = FALSE) 
	        ## Mean number of animals to be tested per holding:
		    nAnimalsMeanPerHerd <- sum(with(nAnimalDataFrame, freq*sampleSize))/length(nAnimalVec)				
	    })
	    nAnimalsMeanPerHerd <- Reduce(function(x,y) c(x,y), nAnimalsMeanPerHerdList)
	    names(nAnimalsMeanPerHerd) <- names(nAnimalsMeanPerHerdList)	
	}	
	return(list(lookupTable = lookupTable, nAnimalsMeanPerHerd = nAnimalsMeanPerHerd))
}
