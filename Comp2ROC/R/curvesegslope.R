curvesegslope <-
function(curve.fpr,curve.tpr){
  curve.segslope=c()
  for (i in 1:(dim(curve.fpr)-1))
  {
    if (curve.fpr[i]!=curve.fpr[i+1])
      curve.segslope[i]=(curve.tpr[i]-curve.tpr[i+1])/(curve.fpr[i]-curve.fpr[i+1])
    else
      curve.segslope[i]=-Inf
  }
  return(curve.segslope)
}
