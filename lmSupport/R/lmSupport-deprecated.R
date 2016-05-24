lm.boxCox <- function (...) {
  warning("'lm.boxCox' is deprecated.\nUse 'modelBoxCox' instead.")
  modelBoxCox(...)
}
lm.codeRegressors <- function (...) {
  warning("'lm.codeRegressors' is deprecated.\nUse 'varRegressors' instead.")
  varRegressors(...)
}

lm.correctSE <- function (...) {
  warning("'lm.correctSE' is deprecated.\nUse 'modelCorrectSE' instead.")
  modelCorrectSE(...)
}

lm.deltaR2 <- function (...) {
	warning("'lm.deltaR2' is deprecated.\nUse 'modelCompare' instead.")
	modelCompare(...)
}

lm.describeGroups <- function (...) {
  warning("'lm.describeGroups' is deprecated.\nUse 'varDescribeBy' instead.")
  varDescribeBy(...)
}

lm.describeData <- function (...) {
  warning("'lm.describeData' is deprecated.\nUse 'varDescribe' instead.")
  varDescribe(...)
}

lm.figSum <- function (...) {
  warning("'lm.figSum' is deprecated.\nUse 'varPlot' instead.")
  varPlot(...)
}

lm.mergeData <- function (...) {
  warning("'lm.mergeData' is deprecated.\nUse 'dfMerge' instead.")
  dfMerge(...)
}

lm.pointEstimates <- function (...) {
  warning("'lm.pointEstimates' is deprecated.\nUse 'modelPredictions' instead.")
  modelPredictions(...)
}  

lm.readDat <- function (...) {
  warning("'lm.readDat' is deprecated.\nUse 'dfReadDat' instead.")
  dfReadDat(...)  
}  

lm.removeCases <- function (...) {
  warning("'lm.removeCases' is deprecated.\nUse 'dfRemoveCases' instead.")
  dfRemoveCases(...)  
} 

lm.renameVar <- function (...) {
  warning("'lm.renameVar' is deprecated.\nUse 'varRename' instead.")
  varRename(...)  
}  

lm.setContrasts <- function (...) {
  warning("'lm.setContrasts' is deprecated.\nUse 'varContrasts' instead.")
  varContrasts(...)  
}   


lm.setRownames <- function (...) {
  warning("'lm.setRownames' is deprecated.\nUse 'dfRownames' instead.")
  dfRownames(...)  
}   

lm.stripChart <- function (...) {
  warning("'lm.stripChart' is deprecated.\nUse 'figStripChart' instead.")
  figStripChart(...)  
}    

lm.sumSquares <- function (...) {
  warning("'lm.sumSquares' is deprecated.\nUse 'modelEffectSizes' instead.")
  modelEffectSizes(...)  
}    

lm.writeDat <- function (...) {
  warning("'lm.writeDat' is deprecated.\nUse 'dfWriteDat' instead.")
  dfWriteDat(...)  
}