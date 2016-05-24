Likelihoods <-
function(resTable){
  sensitivity <- resTable$TP/(resTable$TP + resTable$FN)
  specificity <- resTable$TN/(resTable$FP + resTable$TN)
  
  rhoPos <- sensitivity/(1 - specificity)
  rhoNeg <- (1 - sensitivity)/specificity
  rrho <- cbind(rhoPos = rhoPos, rhoNeg = rhoNeg)
  
  return(rrho)
}
