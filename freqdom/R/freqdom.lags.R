freqdom.lags = function(X){
  if (is.timedom(X))
    lags = X$lags
  else
    lags = X$freq
  lags
}
