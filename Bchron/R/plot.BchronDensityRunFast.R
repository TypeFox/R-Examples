plot.BchronDensityRunFast <-
function(x,plotDates=TRUE,plotSum=FALSE,...) {
  
  # x is the object with everything in it
  # Create a grid over which to plot - the range of the dates plus 10% each side
  n = length(x$calAges)
  nclusters = length(x$clusterRange)
  
  # Find the best cluster method
  chosen = which.max(x$out$BIC)
  
  # Get the means and standard deviations for each chosen cluster
  clusterMeans = x$out$parameters$mean
  clusterSds = sqrt(x$out$parameters$variance$sigmasq)
  clusterProps = x$out$parameters$pro
  
  # Get a range of age values
  thetaRange = round(c(min(clusterMeans-3*clusterSds),max(clusterMeans+3*clusterSds)),0)
  thetaSeq = thetaRange[1]:thetaRange[2]
  
  # Create the densities for each cluster
  dens = matrix(NA,ncol=x$out$G,nrow=length(thetaSeq))
  for(i in 1:ncol(dens)) {
    dens[,i] = stats::dnorm(thetaSeq,mean=clusterMeans[i],sd=clusterSds[i])
  }
    
  # Create a plot
  graphics::plot(thetaSeq,dens%*%clusterProps,type='l',...)
  
  # Add in individual densities
  for(i in 1:ncol(dens)) graphics::lines(thetaSeq,clusterProps[i]*dens[,i],lty=2)
  
  if(plotDates) {
    # Plot the individual dates
    yHeight=graphics::par('usr')[4]
    myCol = grDevices::rgb(190/255,190/255,190/255,0.4)
    for(i in 1:n) {
      graphics::polygon(x$calAges[[i]]$ageGrid,0.3*yHeight*x$calAges[[i]]$densities/max(x$calAges[[i]]$densities),col=myCol,border=NA)
    }   
  }
 
  if(plotSum) {
    yHeight=graphics::par('usr')[4]
    thetaRange = range(x$calAges[[1]]$ageGrid)
    for(i in 2:n) thetaRange = range(c(thetaRange,x$calAges[[i]]$ageGrid))
    dateGrid = seq(round(thetaRange[1]*0.9,0),round(thetaRange[2]*1.1,0),by=1)
    sumDens = rep(0,length(dateGrid))
    for(i in 1:n) {
      matchRows = match(x$calAges[[i]]$ageGrid,dateGrid)
      sumDens[matchRows] = sumDens[matchRows]+x$calAges[[i]]$densities
      if(any(is.na(matchRows))) stop()
    }
    graphics::lines(dateGrid,sumDens*yHeight/max(sumDens),col='red')
  }
  
}
