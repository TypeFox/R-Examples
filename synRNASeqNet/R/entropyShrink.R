entropyShrink <-
function(cellCounts, unit = unit,
                          shrinkageTarget = shrinkageTarget){
  hatTheta <- thetaShrink(cellCounts, shrinkageTarget = shrinkageTarget)
  hatTheta <- hatTheta[hatTheta != 0]
  
  if(unit == "bit") logHatTheta <- log2(hatTheta) else
    if(unit == "ban") logHatTheta <- log10(hatTheta) else
      if(unit == "nat") logHatTheta <- log(hatTheta)
  
  ans <- - sum(hatTheta*logHatTheta)
  return(ans)
}
