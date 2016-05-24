mModeCovariance <-
function(x, m, center = TRUE){
  if(center == TRUE) x <- tensorCentering(x)
  xm <- mFlatten(x, m)
  matrixCovariance(xm)
}
