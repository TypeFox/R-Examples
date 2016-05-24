## Ian Kopacka
## 2010-07-19
##
## Function: computeAlpha
## 
## Compute the herd-based alpha-errors (= 1-herd sensitivity) for a vector
## of herd sizes for limited or individual sampling; see Ziller et al.
##
## Input parameters:
##     nAnimalVec............Numeric vector containing the herd sizes.
##     method................Character string. "individual" for individual sampling
##                           or "limited" for limited sampling.
##     sampleSizeLtd.........Numeric. Sample limit. Required if method == "limited".
##     herdSensitivity.......Numeric between 0 and 1. Desired herd sensitivity. Required
##                           if method == "individual".
##     intraHerdPrevalence...Numeric between 0 and 1. Intra herd prevalence.
##     diagSensitivity.......Numeric between 0 and 1. Sensitivity of the diagnostic
##                           test.
##     diagSpecificity.......Numeric between 0 and 1. Specificity of the diagnostic
##                           test (default value = 1 and is recommended).
##
## Calls: 
##    -
##
## Is called by:
##    method "sample" for classes "IndSampling" and "LtdSampling".
##
computeAlpha <- function(nAnimalVec, method, sampleSizeLtd, herdSensitivity, 
    intraHerdPrevalence, diagSensitivity, diagSpecificity = 1){    
    if (method == "limited"){
        ## Limited sampling:
        ####################
        if (missing(sampleSizeLtd)) stop("Argument 'sampleSizeLtd' must be specified for method == 'limited'")
        
        ## Vector of unique stock sizes:
        nAnimalUniqueVec <- unique(as.numeric(as.character(nAnimalVec)))
        ## Compute corresponding sample sizes (limited sampling):
        sampleSizesUniqueVec <- pmin.int(nAnimalUniqueVec, sampleSizeLtd)         
    } else {
        ## Individual sampling:
        #######################
        if (missing(herdSensitivity)) stop("Argument 'herdSensitivity' must be specified for method == 'individual'")
        ## Vector of unique stock sizes:
        nAnimalUniqueVec <- seq(1,max(nAnimalVec))
        ## Compute corresponding sample sizes (individual sampling):
        sampleSizes <- sapply(nAnimalUniqueVec, 
            function(ii) computeOptimalSampleSize(nPopulation = ii, 
            prevalence = intraHerdPrevalence, alpha = 1-herdSensitivity, 
            sensitivity = diagSensitivity, 
            specificity = diagSpecificity, lookupTable = FALSE))
        nInaccurate <- nAnimalUniqueVec[sampleSizes == Inf]
        sampleSizes[sampleSizes == Inf] <- nInaccurate
        ## Compute cumulative maximum of the sample sizes in order to
        ## force the sample sizes to be monotonically increasing:
        sampleSizesUniqueVec <- rep(1,length(sampleSizes))
        for (ii in seq(along=sampleSizes)[-1]){ 
            sampleSizesUniqueVec[ii] <- max(sampleSizesUniqueVec[ii-1], sampleSizes[ii])
        } 
#		## FOR TESTING ONLY!!!!!!
#		sampleSizesUniqueVec <- sampleSizes
    } 
    ## Compute the number of diseased individuals in the herds (minimum = 1):
    nDiseasedUniqueVec <- round(nAnimalUniqueVec * intraHerdPrevalence)
    nDiseasedUniqueVec <- pmax.int(nDiseasedUniqueVec, 1)
    ## Compute alpha errors for the different herd sizes:
    alphaUniqueVector <- sapply(seq(along = nAnimalUniqueVec), 
        function(ii){computePValue(nPopulation = 
            nAnimalUniqueVec[ii], nSample = sampleSizesUniqueVec[ii], 
            nDiseased = nDiseasedUniqueVec[ii], sensitivity = 
            diagSensitivity, specificity = diagSpecificity)})
    alphaSizesDataFrame <- data.frame(size = nAnimalUniqueVec, 
        alpha = alphaUniqueVector)
    ## Merge with nAnimalVec:
    alphaSizesDataFrame <- merge(x = data.frame(size = nAnimalVec, id = seq(along = nAnimalVec)),
        y = alphaSizesDataFrame, by = "size", sort = FALSE, all = FALSE)
    alphaVec <- alphaSizesDataFrame$alpha[order(alphaSizesDataFrame$id, decreasing = FALSE)]
    ## Return value:
    return(alphaVec)       
}
