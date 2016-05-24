mModeAutoCovariance <-
function(x, m, lag, center=TRUE){
  if(center == TRUE) x <- tensorCentering(x)
  xm <- mFlatten(x, m)
  mAutoCovMatrix(xm, lag)
}
