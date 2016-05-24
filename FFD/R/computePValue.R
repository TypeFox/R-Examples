## Ian Kopacka
## 2010-05-05
##
## Function: computePValue
## 
## Computes the probability of finding no testpositives in
## a sample of a finite population. An imperfect test is used.
##
## Input parameters:
##     nPopulation...Integer. Population size.
##     nSample.......Integer. Sample size.
##     nDiseased.....Integer. Number of diseased in the population (according 
##                   to the null hypothesis).
##     sensitivity...Numeric between 0 and 1. Sensitivity of test (diagnostic 
##                   test for one stage sampling, herd test for two stage 
##                   sampling).
##     specificity...Numeric between 0 and 1. Specificity of test (diagnostic 
##                   test for one stage sampling, herd test for two stage 
##                   sampling).
##
## Source: A.R. Cameron, F.C. Baldock, "A new probability fomula for 
## surveys to substantiate freedom from disease", Prev. Vet. Med. 34
## (1998), pp. 1 - 17.
##
## Calls:
##    -
##
## Is called by:
##    computeOptimalSampleSize.R
##    computeAlphaLimitedSampling.R
##    computeAlpha.R
##
computePValue <- function(nPopulation, nSample, nDiseased, 
        sensitivity, specificity = 1){
    ## Possible number of infected in sample: 
    ##     maximum = min(nSample, nDiseased)
    ## Possible number of healthy in sample: 
    ##     maximum = nPopulation - nDiseased
    ## ==> nSample - (nPopulation - nDiseased) <= nSampleDiseased <= 
    ##          min(nSample, nDiseased)
    nSampleDiseasedVector <- max(0, nSample - (nPopulation - nDiseased)) : 
        min(nSample, nDiseased)
    probabilityHypergoemetricVector <- dhyper(x = nSampleDiseasedVector, 
        m = nDiseased, n = nPopulation - nDiseased, k = nSample)    
    out <- sum(probabilityHypergoemetricVector * 
        (1-sensitivity)^nSampleDiseasedVector * 
        specificity^(nSample-nSampleDiseasedVector))
    return(out)
}
