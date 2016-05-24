DiscriminantPower <-
function(resTable){
  sensitivity <- resTable$TP/(resTable$TP + resTable$FN)
  specificity <- resTable$TN/(resTable$FP + resTable$TN)
  
  X <- sensitivity/(1 - sensitivity)
  Y <- specificity/(1 - specificity)
  DP <- (sqrt(3)/pi) * (log(X) + log(Y))
  
  return(DP)
}
