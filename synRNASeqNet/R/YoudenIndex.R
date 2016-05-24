YoudenIndex <-
function(resTable){
  sensitivity <- resTable$TP/(resTable$TP + resTable$FN)
  specificity <- resTable$TN/(resTable$FP + resTable$TN)
  
  ggamma <- sensitivity - (1 - specificity) #informedness
  return(ggamma)
}
