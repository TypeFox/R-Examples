# Get False Positive Rate(1-Specificity) and True Positive Rate (Sensitivity)
getFprTpr<- function(actuals, predictedScores, threshold=0.5){
  return(list(1-specificity(actuals=actuals, predictedScores=predictedScores, threshold=threshold),
              sensitivity(actuals=actuals, predictedScores=predictedScores, threshold=threshold)))
}
