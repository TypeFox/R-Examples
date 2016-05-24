entropyCS <-
function(cellCounts, unit = unit){
  hatTheta <- thetaGT(cellCounts)
  hatTheta <- hatTheta[hatTheta != 0]
  
  if(unit == "bit") logHatTheta <- log2(hatTheta) else
    if(unit == "ban") logHatTheta <- log10(hatTheta) else
      if(unit == "nat") logHatTheta <- log(hatTheta)
  
  n <- sum(cellCounts)
  
  ans <- - sum((hatTheta*logHatTheta)/(1 - (1 - hatTheta)^n))
  return(ans)
}
