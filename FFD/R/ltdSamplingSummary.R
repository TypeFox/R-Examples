## Ian Kopacka
## 2010-07-14
##
## Function: ltdSamplingSummary
## 
## Constructor for objects of the class 'LtdSamplingSummary'.
## The function computes the mean herd sensitivitiy, the number
## of herds to be tested, the expected total number of animals
## to be tested and the expected total costs for a sequence
## of sample limits (= fixed number of animals to be tested 
## per herd).
##
## Package: FFD
##
## Input parameters:
##    survey.Data........Object of the class 'SurveyData', created by using
##                       the function 'surveyData.R'
##    sampleSizeLtdMax...Positive integer. A series of parameters is computed
##                       for a sequence of sample limits. These sample limits
##                       range from 1 to a given upper bound, defined by
##                       'sampleSizeLtdMax'. If no upper bound is specified
##                       then the maximal herd size is used.
##    nSampleFixVec......Numeric vector containing some NAs (optional argument). 
##                       For risk groups for which the sample size is fixed 
##                       specify the sample size. For the risk groups for which
##                       the sample size should be computed set NA (order of the
##                       risk groups must be the same order as in 
##                       'survey.Data@riskValueData').
##    probVec............Numeric vector. For those risk groups for which the 
##                       sample size should be computed sample probabilities must 
##                       be specified.
##                       The vector must have the same length as the number of 
##                       NA entries in nSampleFixVec or if nSampleFixVec is not 
##                       specified, probVec must have the same length as the 
##                       number of rows in survey.Data@riskValueData.
##
## Return value: object of the class 'LtdSamplingSummary'.

ltdSamplingSummary <- function(survey.Data, sampleSizeLtdMax, 
	nSampleFixVec = NULL, probVec = NULL){
    ## Set value for sampleSizeLtdMax:
	if (missing(sampleSizeLtdMax)) sampleSizeLtdMax <- min(20, max(survey.Data@nAnimalVec))
    if (sampleSizeLtdMax < 1) stop("[ltdSamplingSummary]: sampleSizeLtdMax must be a positive integer.\n")
    if (sampleSizeLtdMax-as.integer(sampleSizeLtdMax) != 0) stop("[ltdSamplingSummary]: sampleSizeLtdMax must be a positive integer.\n")
    sampleSizeLtdVec <- seq(1,sampleSizeLtdMax)    

    outList <- ltdSampling.internal(survey.Data = survey.Data, 
        sampleSizeLtdVec = sampleSizeLtdVec, nSampleFixVec = nSampleFixVec,
		probVec = probVec)
    if (is.null(outList)) return(NULL)
    
    ## Create object of class 'LtdSamplingSummary':
    #if (any(is.na(outList$expectedCostVec))) expectedCostVec <- numeric()
	if (is.null(nSampleFixVec)){
		if (is.null(probVec)){
			nSampleFixVec <- numeric()
		} else {
			nSampleFixVec <- rep(NA, length(probVec))
		    class(nSampleFixVec) <- "numeric"
		}
	} 		
	if (is.null(probVec)) probVec <- numeric()
	
    out <- new("LtdSamplingSummary", 
        surveyData = survey.Data,
        sampleSizeLtdVec = outList$sampleSizeLtdVec,
        meanHerdSensVec = outList$meanHerdSensVec,
		meanHerdSensPerRGMx = outList$meanHerdSensPerRGMx,
        nHerdsVec = outList$nHerdsVec,
		nHerdsPerRiskGroupMx = outList$nHerdsPerRiskGroupMx,
		nSampleFixVec = nSampleFixVec,
		probVec = probVec,
        nAnimalsMeanVec = outList$nAnimalsMeanVec,
        expectedCostVec = outList$expectedCostVec)
    return(out)     
}
