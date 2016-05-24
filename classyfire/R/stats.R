# ************************************************************************
# Main functions for descriptive statistics within classyfire 
# 
# Functions: 
#    getAcc:        Get the test and train accuracies within the ensemble
#    getAvgAcc:     Get the average test and train overall accuracy of the ensemble
#    getOptParam:   Get the optimal hyperparameters
#    getConfMatr:   Get the overall confusion matrix of the ensemble
#    getPerm5Num:   Get the descriptive statistics (five number summary) of the permutation run
#
# ***********************************************************************


# Get the test and train accuracies within the ensemble
getAcc <- function(ensObj) {
  .argsCheck(ensObj, "cfBuild")
  
  testAcc  <- round(ensObj$testAcc,  digits=2)
  trainAcc <- round(ensObj$trainAcc, digits=2)
  
  return (list(Test = testAcc, Train = trainAcc))
}

# Get the average test and train overall accuracy of the ensemble
getAvgAcc  <- function(ensObj) {
  .argsCheck(ensObj, "cfBuild")
  
  avgTest  <- round(mean(ensObj$testAcc),  digits=2)
  avgTrain <- round(mean(ensObj$trainAcc), digits=2)
  
  return (list(Test = avgTest, Train = avgTrain))
}

# Get the optimal hyperparameters 
getOptParam <- function(ensObj) {
  .argsCheck(ensObj, "cfBuild")
  
  gamma     <- ensObj$optGamma
  cost      <- ensObj$optCost
  optHyper  <- cbind(gamma, cost)
  colnames(optHyper) <- c("Opt Gamma", "Opt Cost")
  
  return (optHyper)
}

# Get the overall confusion matrix of the ensemble
getConfMatr <- function(ensObj) {
  .argsCheck(ensObj, "cfBuild")
  
  totalConf <- Reduce("+", ensObj$confMatr)
  propTable <- round(prop.table(totalConf, 1)*100)
  
  return(propTable)
}

# Get the descriptive statistics (five number summary) of the permutation run
getPerm5Num <- function(permObj) {
  .argsCheck(permObj, "cfPermute")
  
  minVal    <- min(permObj$avgAcc)
  maxVal    <- max(permObj$avgAcc)
  medianVal <- median(permObj$avgAcc)
  quantiles <- quantile(permObj$avgAcc)
  
  return (list(minimum = minVal,
               lowerQ  = quantiles[2],
               median  = medianVal,
               upperQ  = quantiles[4],
               maximum = maxVal))
}
