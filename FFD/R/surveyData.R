## Ian Kopacka
## 2010-07-13
##
## Function: surveyData
## 
## Constructor for objects of the class 'SurveyData'.
##
## Package: FFD
##
## Input parameters:
##    nAnimalVec............Integer vector. Stock sizes of the herds.
##    riskGroupVec..........Character vector. Vector containing the 
##                          the name of a risk group to which the farm belongs. 
##                          Optional argument. If provided, it must have the 
#                           same length as nAnimalVec.
##    riskValueData.........Data frame. Data frame where the first column 
##                          contains the labels in riskGroupVec and the second 
##                          column contains the numeric values for the relative 
##                          infection risk.
##    populationData........Data frame. Columns of the data frame
##                          must have the same length as 'nAnimalVec'. 
##                          The data frame can contain additional data such as 
##                          herd id, name and address of the owner etc.
##    designPrevalence......Numeric. Prevalence of the disease under the null hypothesis.
##                          Value between 0 and 1.
##    alpha.................Numeric. Type one error for the statistical test (significance level).
##                          Value between 0 and 1.
##    intraHerdPrevalence...Numeric. Intra-herd prevalence, i.e., the assumed prevalence 
##                          of the disease within an infected herd.
##                          Value between 0 and 1.
##    diagSensitivity.......Numeric. Sensitivity of the diagnostic test.
##                          Value between 0 and 1.
##    costHerd..............Numeric. Cost per tested herd (excluding costs for 
##                          sampling of animals (e.g., travel costs of the vet.))
##    costAnimal............Numeric. Cost per tested animal, e.g., drawing
##                          of samples + analysis in the lab.
##
## Return value: object of the class 'SurveyData'.

surveyData <- function(nAnimalVec = numeric(), riskGroupVec = character(),
	riskValueData = data.frame(),
	populationData = data.frame(), designPrevalence = numeric(), 
	alpha = numeric(), intraHerdPrevalence = numeric(),
    diagSensitivity = numeric(), costHerd = numeric(), costAnimal = numeric()){
    out <- new("SurveyData", nAnimalVec = nAnimalVec, 
		riskGroupVec = riskGroupVec, riskValueData = riskValueData, 
		populationData = populationData,
        designPrevalence = designPrevalence, alpha = alpha, 
        intraHerdPrevalence = intraHerdPrevalence, diagSensitivity = diagSensitivity,
        costHerd = costHerd, costAnimal = costAnimal)
    return(out)    
}
