## ----eval=FALSE----------------------------------------------------------
#  install.packages("classyfire")

## ----eval=FALSE----------------------------------------------------------
#  library(classyfire)

## ----eval=FALSE----------------------------------------------------------
#  ??classyfire

## ----eval=FALSE----------------------------------------------------------
#  data(iris)
#  
#  irisClass <- iris[,5]
#  irisData  <- iris[,-5]

## ----eval=FALSE----------------------------------------------------------
#  ens <- cfBuild(inputData = irisData, inputClass = irisClass, bootNum = 10, ensNum = 10,
#                 parallel = TRUE, cpus = 4, type = "SOCK")

## ----eval=FALSE----------------------------------------------------------
#  ens <- cfBuild(inputData = irisData, inputClass = irisClass, bootNum = 10, ensNum = 10,
#                 parallel = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  attributes(ens)

## ----eval=FALSE----------------------------------------------------------
#  getAvgAcc(ens)$Test
#  getAvgAcc(ens)$Train

## ----eval=FALSE----------------------------------------------------------
#  ens$testAcc
#  ens$trainAcc
#  
#  # Alternatively
#  
#  getAcc(ens)$Test
#  getAcc(ens)$Train

## ----eval=FALSE----------------------------------------------------------
#  testMatr <- matrix(runif(400)*100, ncol = ncol(irisData))
#  predRes  <- cfPredict(ens, testMatr)

## ----eval=FALSE----------------------------------------------------------
#  permObj <- cfPermute(irisData, irisClass, bootNum = 10, ensNum = 10, permNum = 5,
#                       parallel = TRUE, cpus = 4, type = "SOCK")

## ----eval=FALSE----------------------------------------------------------
#  permObj$avgAcc

## ----eval=FALSE----------------------------------------------------------
#  permObj$totalTime[3]
#  permObj$execTime

## ----eval=FALSE----------------------------------------------------------
#  permObj$permList[[1]]

## ----eval=FALSE----------------------------------------------------------
#  getAvgAcc(ens)
#  getAvgAcc(ens)$Test
#  getAvgAcc(ens)$Train

## ----eval=FALSE----------------------------------------------------------
#  getAcc(ens)
#  getAcc(ens)$Test
#  getAcc(ens)$Train

## ----eval=FALSE----------------------------------------------------------
#  getConfMatr(ens)

## ----eval=FALSE----------------------------------------------------------
#  optParam <- getOptParam(ens)
#  optParam

## ----eval=FALSE----------------------------------------------------------
#  getPerm5Num(permObj)
#  getPerm5Num(permObj)$median
#  getPerm5Num(permObj)$minimum
#  getPerm5Num(permObj)$maximum
#  getPerm5Num(permObj)$upperQ
#  getPerm5Num(permObj)$lowerQ

## ----eval=FALSE----------------------------------------------------------
#  # Show the percentages of correctly classified samples in
#  # a barplot with or without text respectively
#  
#  ggClassPred(ens)
#  ggClassPred(ens, showText = TRUE)
#  
#  # Show the percentages of classified and missclassified samples
#  # in a barplot simultaneously with and without text
#  
#  ggClassPred(ens, displayAll = TRUE)
#  ggClassPred(ens, position = "stack", displayAll = TRUE)
#  ggClassPred(ens, position = "stack", displayAll = TRUE, showText = TRUE)
#  
#  # Alernatively, using a dodge position
#  ggClassPred(ens, position = "dodge", displayAll = TRUE)
#  ggClassPred(ens, position = "dodge", displayAll = TRUE, showText = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  ggEnsTrend(ens)
#  
#  # Plot with text
#  ggEnsTrend(ens, showText  = TRUE)
#  
#  # Plot with text; set different limits on y axis
#  ggEnsTrend(ens, showText  = TRUE, ylims=c(90, 100))

## ----eval=FALSE----------------------------------------------------------
#  ggEnsHist(ens)
#  
#  # Density plot of the test accuracies in the ensemble
#  ggEnsHist(ens, density = TRUE)
#  
#  # Density plot that highlights additional descriptive statistics
#  ggEnsHist(ens, density = TRUE, percentiles=TRUE)
#  ggEnsHist(ens, density = TRUE, percentiles=TRUE, mean=TRUE)
#  ggEnsHist(ens, density = TRUE, percentiles=TRUE, median=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  ggPermHist(permObj)
#  
#  # Density plot
#  ggPermHist(permObj, density=TRUE)
#  
#  # Density plot that highlights additional descriptive statistics
#  ggPermHist(permObj, density=TRUE, percentiles = TRUE, mean = TRUE)
#  ggPermHist(permObj, density=TRUE, percentiles = TRUE, median = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  ggFusedHist(ensObj, permObj)

