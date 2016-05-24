plot.RescorlaWagner = function(x, asymptote=TRUE, xlab="t", ylab="weight",
    ylimit=NA, ...)
{
  rwWeights = x 
  x = 1:length(rwWeights$weightvector)
  y = rwWeights$weightvector
  eqw = rwWeights$equilibriumWeight
  if (is.na(ylimit[1])) {
    if (asymptote) {
      ylimit = range(y, eqw)
    } else {
      ylimit = range(y)
    }
  } 
  plot.default(x, y, type="l", ylim=ylimit, xlab=xlab, ylab=ylab, ...)
  if (asymptote) abline(h=eqw, lty=2, ...)
}
