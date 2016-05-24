# **************************************************************************************************************
# Functions for unit testing
# **************************************************************************************************************

library('RUnit')
library('classyfire')

set.seed(1)

# Test data
data(iris)
irisClass <- iris[,5]
irisData  <- iris[,-5]
randClass <- c(2, rep(3, length(irisClass)-1))
testVec  <- t(c(62,20,68,76))      

# Use parallel = FALSE for testing on CRAN! 
ensObj  <- cfBuild(inputData = irisData, inputClass = irisClass, bootNum = 5, ensNum = 2, parallel = FALSE)
permObj <- cfPermute(irisData, irisClass, bootNum = 5, ensNum = 2, permNum = 2, parallel = FALSE)
predRes <- cfPredict(ensObj , testVec)


# Test the initial checks on the input data provided by the user etc. 
test.initCheck <- function() {
  checkException(.initCheck(), silent=TRUE)
  checkException(.initCheck(irisData), silent=TRUE)
  checkException(.initCheck(inputClass = inputClass), silent=TRUE)
  checkException(.initCheck(iris, irisClass),     silent=TRUE)
  checkException(.initCheck(irisData, irisData),  silent=TRUE)
  checkException(.initCheck(irisData, irisData),  silent=TRUE)
  checkException(.initCheck(irisData, randClass), silent=TRUE)
}

# Test the main cfBuild function for the construction of the ensemble
test.cfBuild <- function() {
  checkEquals("cfBuild",  class(ensObj)[2])
  checkEquals(13,         length(ensObj))
  checkEquals(95.1,       getAvgAcc(ensObj)$Test)
  checkEquals(96.97,      getAvgAcc(ensObj)$Train) 
  checkEquals(92.16,      ensObj$testAcc[1])
  checkEquals(96.97,      ensObj$trainAcc[1])
  checkEquals(TRUE,       any(attributes(ensObj)$names == "testAcc"))
  checkEquals(100,        getConfMatr(ensObj)[1,1])
  checkEquals(94,         getConfMatr(ensObj)[2,2])
  checkEquals(91,         getConfMatr(ensObj)[3,3])
}

# Test the cfPermute function for permutation testing
test.cfPermute <- function() {
  checkEquals("cfPermute",  class(permObj)[2])
  checkEquals(4,            length(permObj))
  checkEquals(39.22,        permObj$avgAcc[1])
  checkEquals(2,            length(permObj$permList))
}

# Test the cfPredict function for use with unknown data
test.cfPredict <- function() {
  checkEquals("virginica",  as.character(predRes[,1]))
  checkEquals(100,          predRes[,2])
}

# Test the relevant stats functions
test.stats <- function() { 
  checkException(getAcc(),      silent=TRUE)
  checkException(getAvgAcc(),   silent=TRUE)
  checkException(getOptParam(), silent=TRUE)
  checkException(getConfMatr(), silent=TRUE)
  checkException(getPerm5Num(), silent=TRUE)
  
  checkException(getAcc(randClass),      silent=TRUE)
  checkException(getAvgAcc(randClass),   silent=TRUE)
  checkException(getOptParam(randClass), silent=TRUE)
  checkException(getConfMatr(randClass), silent=TRUE)
  checkException(getPerm5Num(randClass), silent=TRUE)
  
  checkEquals(2,         length(getAcc(ensObj)))
  checkEquals(2,         length(getAvgAcc(ensObj)))
  checkEquals(92.16,     getAcc(ensObj)$Test[1])
  checkEquals(96.97,     getAcc(ensObj)$Train[1])
  checkEquals(95.10,     getAvgAcc(ensObj)$Test)
  checkEquals(96.97,     getAvgAcc(ensObj)$Train)
  checkEquals("matrix",  class(getOptParam(ensObj)))
  checkEquals("table",   class(getConfMatr(ensObj)))
  checkEquals(9,         length(getConfMatr(ensObj)))
  checkEquals(5,         length(getPerm5Num(permObj)))
  checkEquals(33.33,     getPerm5Num(permObj)$minimum)
}

# Test the relevant plot functions
test.plots <- function() { 
  checkException(ggEnsTrend(),   silent=TRUE)
  checkException(ggEnsHist(),    silent=TRUE)
  checkException(ggClassPred(),  silent=TRUE)
  checkException(ggPermHist(),   silent=TRUE)
  
  checkException(ggEnsTrend(permObj),   silent=TRUE)
  checkException(ggEnsHist(permObj),    silent=TRUE)
  checkException(ggClassPred(permObj),  silent=TRUE)
  checkException(ggPermHist(ensObj),    silent=TRUE)
}

# Execute all the tests
test.initCheck()
test.cfBuild()
test.cfPredict()
test.cfPermute()
test.stats()
test.plots()
