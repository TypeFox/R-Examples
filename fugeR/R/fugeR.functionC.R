#mfId : the membership function to use
#mfPoints : starting point of each membership function
#values : data to fuzzify
#minValues : As we apply only the operator MIN (AND) we only
#            keep the min values as rule activation.
cfuzzifyVar <- function( mfId, mfPoints, values, minValues ) {
    .Call('cfuzzifyVar', mfId, mfPoints, values, minValues, PACKAGE='fugeR')
}


#nbCase : Number of case in the training data
#nbRule : Number of rule in the fuzzySystem
#inRule : Logical values indicating which rules use the ouput var
#lstMf  : The list of the singleton for this out var
#lstMfId : The singleton used by each rule
#defaultMfId : The singleton used in the default rule
#lstActivation : The activation of each rule for each case in the data.
cdefuzzify <- function( nbCase, nbRule, inRule, lstMf,
                        lstMfId, defaultMfId, lstActivation ) {
    .Call('cdefuzzify', nbCase, nbRule, inRule, lstMf,
                        lstMfId, defaultMfId, lstActivation, PACKAGE='fugeR')
}

#lstPredicted : Values predicted by the model
#lstActual : Actual values
#threshold : threshold to apply on actual and predicted values.
cbinaryPerformance <- function( lstPredicted, lstActual, threshold ) {
    .Call('cbinaryPerformance', lstPredicted, lstActual, threshold, PACKAGE='fugeR')
}


