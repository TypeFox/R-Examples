thetaShrink <-
function(cellCounts, shrinkageTarget = shrinkageTarget){
  hatThetaML <- thetaML(cellCounts)
  n <- sum(cellCounts)
  
  if(is.null(shrinkageTarget)){
    cellCounts <- as.matrix(cellCounts)
    p <- nrow(cellCounts)*ncol(cellCounts)
    shrinkageTarget <- 1/p
  }
  
  lambda <- shrinkageIntensity(hatThetaML = hatThetaML, n = n,
                               shrinkageTarget = shrinkageTarget)
  
  ans <- lambda*shrinkageTarget + (1 - lambda)*hatThetaML
  return(ans)
}
