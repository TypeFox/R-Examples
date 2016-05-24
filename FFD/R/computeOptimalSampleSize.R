## Ian Kopacka
## 2010-05-05
##
## Function: computeOptimalSampleSize
## 
## Computes the optimal sample size for a survey to substantiate
## freedom from disease. The optimal sample size is the smallest
## sample size that produces an alpha-error less than or equal 
## to a prediscribed value for alpha. The population is considered 
## as diseased if at least one individual has a positive test result.
## The sample size is computed using a bisection method. 
## The function either computes the sample size for a fixed population
## (lookupTable = FALSE) or a lookup table with sampling sizes 
## depending on the population size for individual sampling
## (lookupTable = TRUE); see Ziller et al.
##
## Input parameters:
##     nPopulation...Integer. Population size.
##     prevalence....Numeric between 0 and 1. Prevalence of the disease 
##                   under the null hypothesis.
##     alpha.........Numeric between 0 and 1. Type one error for the 
##                   statistical test.
##     sensitivity...Numeric between 0 and 1. Sensitivity of test (diagnostic 
##                   test for one stage sampling, herd test for two stage 
##                   sampling).
##     specificity...Numeric between 0 and 1. Specificity of test (diagnostic 
##                   test for one stage sampling, herd test for two stage 
##                   sampling).
##     lookupTable...Logical. TRUE if a lookup table of sample sizes for
##                   individual sampling (see, Ziller et al., 2002) should
##                   be produced. FALSE if the sample size is desired
##                   for a fixed population size (default).
##
## Calls: 
##    computePValue.R
##
## Is called by:
##    -
##
computeOptimalSampleSize <- function(nPopulation, prevalence, alpha = 0.05, 
    sensitivity = 1, specificity = 1, lookupTable = FALSE){
    if(!lookupTable){
        ## Compute sample size for one fixed population size:
        #####################################################
           
        ## Calculate the number of diseased individuals in the population:
        nDiseased <- max(round(nPopulation*prevalence),1)
        ## Initialize the parameters for the bisection method:
        lowerBoundSampleSize <- 1
        upperBoundSampleSize <- nPopulation
        probabilityLowerUpper <- sapply(c(lowerBoundSampleSize, 
            upperBoundSampleSize), 
            function(sampleSize) computePValue(nPopulation = 
            nPopulation, nSample = sampleSize, nDiseased = nDiseased, 
            sensitivity = sensitivity, specificity = specificity))
        ## Check if upper bound satisfies the desired accuracy:
        if (probabilityLowerUpper[2] > alpha) 
            return(Inf)
        ## Check if lower bound satisfies the desired accuracy:
        if (probabilityLowerUpper[1] <= alpha) 
            return(lowerBoundSampleSize)
        
        ## Bisection method: Due to the discrete nature of the variable
        ## (= optimal sample size) the bisection method is executed
        ## until the width of the search interval falls below a certain
        ## threshold. The probabilities for the remaining cases are 
        ## explicitly computed and the optimal sample size is chosen:
        minIntervalWidth <- 4 #10
        while (upperBoundSampleSize - lowerBoundSampleSize > minIntervalWidth){
            midSampleSize <- round((lowerBoundSampleSize + upperBoundSampleSize)/2)
            probabilityMid <- computePValue(nPopulation = 
                nPopulation, nSample = midSampleSize, nDiseased = nDiseased, 
                sensitivity = sensitivity, specificity = specificity)
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
            function(sampleSize) computePValue(nPopulation = 
            nPopulation, nSample = sampleSize, nDiseased = nDiseased, 
            sensitivity = sensitivity, specificity = specificity))
        indexProbablityAlpha <- which(probabilityVector <= alpha)
        if (length(indexProbablityAlpha) > 0){
            out <- sampleSizeVector[min(indexProbablityAlpha)]
        } else {
            out <- upperBoundSampleSize
        }
        return(out)
    } else {
        ## Create lookup table for individual sampling, i.e., calculate
        ## sampling sizes for a series of population sizes and return
        ## sample sizes in the form of a lookup table.
        populationSizes <- seq(1,nPopulation)
        sampleSizes <- sapply(populationSizes, 
            function(ii) computeOptimalSampleSize(nPopulation = ii, 
            prevalence = prevalence, alpha = alpha, sensitivity = sensitivity, 
            specificity = specificity, lookupTable = FALSE))
        nInaccurate <- populationSizes[sampleSizes == Inf]
        if(length(nInaccurate) > 0){
            sampleSizes[sampleSizes == Inf] <- nInaccurate            
            warning(paste("Desired herd sensitivity could not be achieved for population sizes ", 
                range(nInaccurate)[1], " - ", range(nInaccurate)[2], ".", sep = ""))
        }
        ## Compute cumulative maximum of the sample sizes in order to
        ## force the sample sizes to be monotonically increasing:
        cumMaxSampleSizes <- rep(1,nPopulation)
        for (ii in 2:nPopulation) cumMaxSampleSizes[ii] <- max(cumMaxSampleSizes[ii-1], sampleSizes[ii])
        ## Create lookup table:
        splitPopulation <- split(x = populationSizes, f = cumMaxSampleSizes)
        rangePopulation <- lapply(splitPopulation, function(x) range(x))
        out <- Reduce(function(x,y) rbind(x,y), rangePopulation)
        if (!is.matrix(out)){
            out <- matrix(out,nrow = 1)
        }        
        out <- cbind(out, as.numeric(names(rangePopulation)))
        colnames(out) <- c("N_lower", "N_upper", "sampleSize")
        rownames(out) <- NULL
        return(out)    
    }
}
