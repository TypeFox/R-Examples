curvesegsloperef <-
function(curve.fpr,curve.tpr,ref.point) {
  curve.slope=c()
  for (i in 1:dim(curve.fpr))
  {
    if (curve.fpr[i]-ref.point[1]==0)
      curve.slope[i]=-Inf
    else
      curve.slope[i]=(curve.tpr[i]-ref.point[2])/(curve.fpr[i]-ref.point[1])
  }
  return(curve.slope)
}
