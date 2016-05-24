summary.value <-
function (object, ...){
  
  res <- list (optVal=object$value, PosPos=object$valPosPos,
  PosNeg=object$valPosNeg, NegPos=object$valNegPos,
  NegNeg=object$valNegNeg); 
  class (res) <- "summary.value"
  res
}
