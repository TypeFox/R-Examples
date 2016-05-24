## functions to calculate performance

calcPerformance = function(annotations, predictions, winSize, names=NULL, 
                           combineStanding=FALSE) {
  cat("\n")
  labelDir = annotationsToLabels(annotations, winSize, names)
  calcPerformanceFromLabels(labelDir, predictions, names, combineStanding)
}
calcPerformanceFromLabels = function(labelDir, predDir, names=NULL, 
                                     combineStanding=FALSE) {
  data = loadPredictionsAndLabels(labelDir, predDir, names)
  if (nrow(data) > 0) {
    pr = data[data$behavior != "NULL", c("prediction")]
    gt = data[data$behavior != "NULL", c("behavior")]
    if (combineStanding) {
      pr[grepl("Standing", pr)] = "Standing"
      gt[grepl("Standing", gt)] = "Standing"
    }
    l = unique(c(pr, gt))
    pr = factor(pr, levels=l)
    gt = factor(gt, levels=l)
    m = confusionMatrix(pr, gt)
    print(m)
    cat("\n")
    return(m)
  }
  else {
    return(NULL)
  }
}