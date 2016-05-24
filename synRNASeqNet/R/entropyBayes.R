entropyBayes <-
function(cellCounts, unit = unit,
                         priorHyperParam = priorHyperParam){
  hatTheta <- thetaBayes(cellCounts, priorHyperParam = priorHyperParam)
  p0 <- hatTheta$p - prod(dim(hatTheta$thetak))
  hatTheta <- c(hatTheta$thetak, hatTheta$theta0)
  
  if(unit == "bit") logHatTheta <- log2(hatTheta) else
    if(unit == "ban") logHatTheta <- log10(hatTheta) else
      if(unit == "nat") logHatTheta <- log(hatTheta)
  hatTheta[length(hatTheta)] <- p0*hatTheta[length(hatTheta)]
  
  ans <- - sum(hatTheta*logHatTheta)
  return(ans)
}
