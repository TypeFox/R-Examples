## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library("knitr")
render_sweave()
# set global chunk options
opts_chunk$set(prompt=TRUE)
options(replace.assign=TRUE, width=90, prompt="R> ")

## ----echo=FALSE, results='hide', message=FALSE------------------------------------------
library("RJafroc")

## ---------------------------------------------------------------------------------------
str(rocData)

## ---------------------------------------------------------------------------------------
str(frocData)

## ---------------------------------------------------------------------------------------
str(roiData)

## ---------------------------------------------------------------------------------------
rocXlsx <- system.file("tests", "rocData.xlsx", package = "RJafroc")
rocLrc <- system.file("tests", "rocData.lrc", package = "RJafroc")
rocCsv <- system.file("tests", "rocData.csv", package = "RJafroc")
rocImrmc <- system.file("tests", "rocData.imrmc", package = "RJafroc")
frocXlsx <- system.file("tests", "frocData.xlsx", package = "RJafroc")
roiXlsx <- system.file("tests", "roiData.xlsx", package = "RJafroc")

RocDataXlsx<- ReadDataFile(fileName = rocXlsx)
RocDataLrc<- ReadDataFile(fileName = rocLrc, format = "MRMC")
RocDataCsv<- ReadDataFile(fileName = rocCsv, format = "MRMC")
RocDataImrmc<- ReadDataFile(fileName = rocImrmc, format = "iMRMC")
FrocDataXlsx<- ReadDataFile(fileName = frocXlsx)
RoiDataXlsx<- ReadDataFile(fileName = roiXlsx)


## ----results='hide'---------------------------------------------------------------------
retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
print(retDbmRoc)


## ----results='hide'---------------------------------------------------------------------
retORRoc  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
print(retORRoc)
CovOR <- retORRoc$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varEps <- CovOR$varCov[6]
msTR <- retORRoc$msTR
msT <- retORRoc$msT

## ---------------------------------------------------------------------------------------
print(CovOR)

## ---------------------------------------------------------------------------------------
retDbm  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
varYTR <- retDbm$varComp$varComp[3]
varYTC <- retDbm$varComp$varComp[4]
varYEps <- retDbm$varComp$varComp[6]

## ---------------------------------------------------------------------------------------
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, 
    effectSize = effectSize) 
  message("# of readers = ", J, ", estimated # of cases = ", ret$K, "\n",
    "predicted power = ", signif(ret$power, 4), "\n")
}

## ---------------------------------------------------------------------------------------
retOR  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
CovOR <- retOR$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varErrOR <- CovOR$varCov[6]
msTR <- retOR$msTR
KStar <- length(rocData$NL[1,1,,1])
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, cov1 = cov1, cov2 = cov2, cov3 = cov3, 
    varEps = varErrOR, msTR = msTR, KStar = KStar, effectSize = effectSize) 
  message("# of readers = ", J, ", estimated # of cases = ", ret$K, "\n",
    "predicted power = ", signif(ret$power, 4), "\n")
}

## ----results='hide',eval=FALSE----------------------------------------------------------
#  retDbmwJafroc  <- DBMHAnalysis(frocData)
#  print(retDbmwJafroc)

## ----results='hide', eval=FALSE---------------------------------------------------------
#  retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")
#  print(retDbmwJafroc1)
#  
#  retDbmJafroc  <- DBMHAnalysis(frocData, fom = "JAFROC")
#  print(retDbmJafroc)
#  
#  retDbmJafroc1  <- DBMHAnalysis(frocData, fom = "JAFROC1")
#  print(retDbmJafroc1)

## ----results='hide',eval=FALSE----------------------------------------------------------
#  retDbmHrAuc  <- DBMHAnalysis(frocData, fom = "HrAuc")
#  retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1")
#  retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2")

## ----results='hide',eval=FALSE----------------------------------------------------------
#  retDbmRoi  <- DBMHAnalysis(roiData, fom = "ROI")

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = rocData, fom = "Wilcoxon",
#    dataDscrpt = "MyROCData", showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = rocData, fom = "Wilcoxon",
#    reportFile = "MyROCDataAnalysis.txt",  showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = rocData, method = "ORH", fom = "Wilcoxon",
#    showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = frocData, fom = "Wilcoxon",
#    showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = frocData, method = "ORH",
#    showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = frocData, fom = "HrAuc",
#    showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport(dataset = roiData, method = "ORH", fom = "ROI",
#    showWarnings = FALSE)

## ----message=FALSE,eval=FALSE-----------------------------------------------------------
#  OutputReport("rocData.xlsx", format = "JAFROC", method = "DBMH",
#    fom = "Wilcoxon", dataDscrpt = "MyROC2Data",  showWarnings = FALSE)

## ---------------------------------------------------------------------------------------
SaveDataFile(dataset = rocData, fileName = "rocData2.xlsx", 
  format = "JAFROC")
SaveDataFile(dataset = rocData, fileName = "rocData2.csv", format = "MRMC")
SaveDataFile(dataset = rocData, fileName = "rocData2.lrc", format = "MRMC", 
  dataDscrpt = "ExampleROCdata")
SaveDataFile(dataset = rocData, fileName = "rocData2.txt", format = "MRMC", 
  dataDscrpt = "ExampleROCdata2")
SaveDataFile(dataset = rocData, fileName = "rocData.imrmc", 
  format = "iMRMC", dataDscrpt = "ExampleROCdata3") 

## ---------------------------------------------------------------------------------------
plotM <- c(1:2)
plotR <- c(1:5)
plotROC <- EmpiricalOpCharac(data = rocData, trts = plotM, rdrs = plotR, 
  opChType = "ROC")

## ---------------------------------------------------------------------------------------
plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:5),c(1:5))
plotRocAvg <- EmpiricalOpCharac(dataset = rocData, trts = plotMAvg, 
  rdrs = plotRAvg, opChType = "ROC")

## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotROC$ROCPlot)

## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotRocAvg$ROCPlot)

## ---------------------------------------------------------------------------------------
plotM <- c(1:2)
plotR <- c(1:4)
plotROC <- EmpiricalOpCharac(data = frocData, trts = plotM, rdrs = plotR, 
  opChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotRocAvg <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
  rdrs = plotRAvg, opChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotAFROC <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
  rdrs = plotRAvg, opChType = "AFROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotFROC <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
  rdrs = plotRAvg, opChType = "FROC")


## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotROC$ROCPlot)

## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotRocAvg$ROCPlot)

## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotAFROC$AFROCPlot)

## ----out.width= '0.49\\linewidth', fig.show='hold', echo=FALSE--------------------------
print(plotFROC$FROCPlot)

## ----eval=FALSE-------------------------------------------------------------------------
#  RJafrocGui()

## ----eval=FALSE-------------------------------------------------------------------------
#  RJafrocGui(useBrowser = TRUE)

