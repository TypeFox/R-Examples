entropyML <-
function(cellCounts, unit = unit){
  hatTheta <- thetaML(cellCounts)
  hatTheta <- hatTheta[hatTheta != 0]
  
  if(unit == "bit") logHatTheta <- log2(hatTheta) else
    if(unit == "ban") logHatTheta <- log10(hatTheta) else
      if(unit == "nat") logHatTheta <- log(hatTheta)# else
  #stop("Unknown Entropy Unit") #match.arg
  
  ans <- - sum(hatTheta*logHatTheta)
  return(ans)
}
