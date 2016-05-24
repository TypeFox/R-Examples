## ----setup, include=FALSE---------------------------------
library(rmarkdown)
options(continue=" ")
options(width=60)
library(knitr)
library(EGRET)
library(EGRETci)


## ----eval=FALSE, echo=TRUE--------------------------------
#  library(EGRET)
#  library(EGRETci)
#  eList <- Choptank_eList
#  
#  #Interactive function to set up trend analysis:
#  caseSetUp <- trendSetUp(eList)
#  eBoot <- wBT(eList,caseSetUp,
#               fileName = "outputText.txt")
#  
#  #Interactive save output function:
#  saveEGRETci(eList, eBoot, caseSetUp)
#  
#  

## ----eval=FALSE-------------------------------------------
#  library(EGRET)
#  library(EGRETci)
#  eList <- Choptank_eList
#  eList <- setPA(eList, paStart = 12, paLong = 4)
#  eList$INFO$windowY <- 10
#  eList$INFO$minNumObs <- 50
#  caseSetUp <- trendSetUp(eList,
#                          year1=1990,
#                          year2=2012,
#                          nBoot = 50,
#                          bootBreak = 39,
#                          blockLength = 200)
#  eBoot <- wBT(eList, caseSetUp, fileName ="outputText.txt")
#  saveEGRETci(eList, eBoot, caseSetUp, fileName = "output")

## ---- fig.height=6, fig.width=6---------------------------
library(EGRET)
library(EGRETci)

# Example data included in package:
eList <- Choptank_eList # Example data from EGRET package
eBoot <- Choptank_eBoot
caseSetUp <- Choptank_caseSetUp

#Concentration an initial run:
plotHistogramTrend(eList, eBoot, caseSetUp,  
                   flux=FALSE)

#Flux an initial run:
plotHistogramTrend(eList, eBoot, caseSetUp,
                   flux=TRUE)


## ---- histExampleCombo, fig.width=7, fig.height=4---------
par(mfrow=c(1,2))
#Concentration, presentation version:
plotHistogramTrend(eList, eBoot, caseSetUp,  
                   flux=FALSE, xMin = -5, xMax = 65, xStep = 5)

#Flux, presentation version:
plotHistogramTrend(eList, eBoot, caseSetUp,
                   flux=TRUE, xMin = -5, xMax = 55, xStep = 5)
# To return to figures printing in 1 row, 1 columns:
par(mfrow=c(1,1))                   

## ---- histExampleCombo2, fig.width=7, fig.height=10-------
par(mfrow=c(2,1))
#Concentration, presentation version:
plotHistogramTrend(eList, eBoot, caseSetUp,  
                   flux=FALSE, xMin = -5, xMax = 65, xStep = 5)

#Flux, presentation version:
plotHistogramTrend(eList, eBoot, caseSetUp,
                  flux=TRUE, xMin = -5, xMax = 65, xStep = 5,
                  printTitle = FALSE)
# To return to figures printing in 1 row, 1 columns:
par(mfrow=c(1,1))                   

## ---- eval=FALSE, echo=TRUE-------------------------------
#  eList <- modelEstimation(eList)

## ---- eval=FALSE------------------------------------------
#  library(EGRET)
#  library(EGRETci)
#  
#  eList <- Choptank_eList
#  
#  CIAnnualResults <- ciCalculations(eList)
#  
#  save(eList,CIAnnualResults, file="CIAnnualResults.RData")

## ----eval=FALSE-------------------------------------------
#  CIAnnualResults <- ciCalculations(eList, nBoot = 100, blockLength = 200, widthCI = 90)

## ---- eval=FALSE------------------------------------------
#  library(foreach)
#  library(doParallel)
#  library(iterators)
#  library(EGRET)
#  library(EGRETci)
#  
#  eList <- Choptank_eList
#  eList <- modelEstimation(eList)
#  
#  nBoot <- 100
#  blockLength <- 200
#  coreOut <- 1 #Number of cores to leave out of processing tasks
#  
#  widthCI <- 90
#  ciLower <- (50-(widthCI/2))/100
#  ciUpper <- (50+(widthCI/2))/100
#  probs <- c(ciLower,ciUpper)
#  
#  nCores <- detectCores() - coreOut
#  cl <- makeCluster(nCores)
#  registerDoParallel(cl)
#  repAnnual <- foreach(n = 1:nBoot,.packages=c('EGRETci')) %dopar% {
#     annualResults <- bootAnnual(eList, blockLength)
#  }
#  stopCluster(cl)
#  
#  # save(repAnnual, file="repAnnual.RData")
#  
#  CIAnnualResults <- ciBands(eList, repAnnual, probs)
#  save(eList,CIAnnualResults, file="CIAnnualResults.RData")
#  

## ---- fig.height=5, fig.width=7---------------------------
eList <- Choptank_eList

CIAnnualResults <- Choptank_CIAnnualResults

plotConcHistBoot(eList, CIAnnualResults)

plotFluxHistBoot(eList, CIAnnualResults)


