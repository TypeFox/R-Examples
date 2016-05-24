plot.BchronologyRun <-
function(x,...) {
  # x contains the output from a run of the Bchronology function
  
  # Get chronology ranges
  chronLow = apply(x$thetaPredict,2,'quantile',probs=0.025)
  chronMed = apply(x$thetaPredict,2,'quantile',probs=0.5)
  chronHigh = apply(x$thetaPredict,2,'quantile',probs=0.975)
  
  xLimits = range(c(chronLow,chronHigh))
  yLimits = range(x$predictPositions)
  for(i in 1:length(x$calAges)) {
    yLimits = range(c(yLimits,x$calAges[[i]]$positions))
  }
  yLimits = c(yLimits[1],yLimits[2])
  dateHeight=0.2*diff(pretty(yLimits))[1]
  yLimits[1] = yLimits[1]-dateHeight
  
  # Create the plot
  graphics::plot(1,1,xlim=rev(xLimits),ylim=rev(yLimits),type="n",...)
  graphics::grid()
  
  # Add in the dates
  for(i in 1:length(x$calAges)) {
    # First for known points
    if(x$calAges[[i]]$ageSds<5) {
      graphics::points(sum(x$calAges[[i]]$ageGrid*x$calAges[[i]]$densities),x$calAges[[i]]$positions,pch=16)
    } else {
      graphics::polygon(c(min(x$calAges[[i]]$ageGrid),x$calAges[[i]]$ageGrid,max(x$calAges[[i]]$ageGrid)),c(x$calAges[[i]]$positions,x$calAges[[i]]$positions-x$calAges[[i]]$dens*dateHeight/max(x$calAges[[i]]$dens),x$calAges[[i]]$positions),border=NA,col='black')      
    }
  }
  
  # Add in the chronologies
  chronCol = grDevices::rgb(190/255,190/255,190/255,alpha=0.8)
  chronLow = apply(x$thetaPredict,2,'quantile',probs=0.025)
  chronMed = apply(x$thetaPredict,2,'quantile',probs=0.5)
  chronHigh = apply(x$thetaPredict,2,'quantile',probs=0.975)
  graphics::polygon(c(chronLow,rev(chronHigh)),c(x$predictPositions,rev(x$predictPositions)),col=chronCol,border=NA)

  graphics::legend('topleft',c('Dated positions','95% Chronology CI'),col=c('black',chronCol),pch=15)
  
}
