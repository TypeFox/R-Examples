## Ian Kopacka
## 2010-07-16
##
## Function: indSamplingSummary
## 
## Constructor for objects of the class 'IndSamplingSummary'.
## The function computes the number of herds to be tested, 
## the expected total number of animals to be tested and the 
## expected total costs for a sequence of herd sensitivities
## ranging from 0.1 to the sensitivity of the diagnostic test.
## The step size for the herd sensitivities can be specified
## by the user. If no step size is specified a step size of
## 0.02 is used.
##
## Package: FFD
##
## Input parameters:
##    survey.Data.....Object of the class 'SurveyData', created by using
##                    the function 'surveyData.R'
##    stepSize........Numeric. A series of parameters is computed for a 
##                    sequence of herd sensitivities. The argument 'stepSize'
##                    specifies the step size used in the discretization of
##                    the herd sensitivities (default = 0.02).
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
## Return value: object of the class 'IndSamplingSummary'.

indSamplingSummary <- function(survey.Data, stepSize = 0.02, 
	nSampleFixVec = NULL, probVec = NULL){
    ## Set value for stepSize:
    herdSensMin <- 0.1
    if (stepSize > (survey.Data@diagSensitivity - herdSensMin)) stop("[indSamplingSummary]: stepSize exceeds range of herd sensitivites.\n")
    if (stepSize <= 0) stop("[indSamplingSummary]: stepSize must be positive.\n")
    herdSensVec <- seq(herdSensMin, survey.Data@diagSensitivity, by = stepSize)    

    outList <- indSampling.internal(survey.Data = survey.Data, 
        herdSensVec = herdSensVec, nSampleFixVec = nSampleFixVec,
		probVec = probVec)    
    if (is.null(outList)) return(NULL)
	
	if (is.null(nSampleFixVec)){
		if (is.null(probVec)){
			nSampleFixVec <- numeric()
		} else {
			nSampleFixVec <- rep(NA, length(probVec))
		    class(nSampleFixVec) <- "numeric"
		}
	} 		
	if (is.null(probVec)) probVec <- numeric()

    ## Create object of class 'IndSamplingSummary':
    if (any(is.na(outList$expectedCostVec))) expectedCostVec <- numeric()
    out <- new("IndSamplingSummary", 
        surveyData = survey.Data,
        herdSensVec = outList$herdSensVec,        
        nHerdsVec = outList$nHerdsVec,
		nHerdsPerRiskGroupMx = outList$nHerdsPerRiskGroupMx,
		nSampleFixVec = nSampleFixVec,
		probVec = probVec,
        nAnimalsMeanVec = outList$nAnimalsMeanVec,
        expectedCostVec = outList$expectedCostVec)
    return(out)     
}
