plotDataPPC <-
function(y, mu, sigma, nu, xVec, maxY) {
  # Does the plots of posterior predictive curves for one sample; called by plotAll;
  #  not exported.
  # Does not do title or sample size: those are added later.
  # y = original data for this sample; can be NULL
  # mu, sigma, nu = vectors of parameters to use for the t-curves and histogram breaks.
  # xVec = vector of values to use for the x-axis values
  # maxY = height of the y-axis

  plot(xVec[1], 0, xlim=range(xVec), ylim=c(0, maxY), cex.lab=1.75,
        type="n", xlab="y", ylab="p(y)")
  for ( i in seq_along(mu)) {
    lines(xVec, dt( (xVec-mu[i])/sigma[i], df=nu[i] )/sigma[i], col="skyblue")
  }

  if(!is.null(y)) {
    histBinWd <- median(sigma)/2
    histCenter <- mean(mu)
    histBreaks <- try(sort( c( seq( histCenter-histBinWd/2 , min(xVec)-histBinWd,#/2 ,
                             -histBinWd ),
                        seq( histCenter+histBinWd/2 , max(xVec)+histBinWd,#/2 ,
                             histBinWd ) ) ), silent=TRUE)
    if(inherits(histBreaks, "try-error")) {
      text(min(xVec), maxY, "Cannot plot data.", pos=4, col='red')
    } else {
      histInfo <- hist( y, plot=FALSE , breaks=histBreaks )
      PlotMat <- cbind(histInfo$mids, histInfo$density)
      PlotMat[histInfo$density == 0] <- NA
      points( PlotMat, type="h" , lwd=3 , col="red" )
    }
  }
}
