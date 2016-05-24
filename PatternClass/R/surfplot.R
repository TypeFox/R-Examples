surfplot <-
function(metric=9, prop=0.7, rho=0.2, colour=TRUE, drop=TRUE, cross=TRUE, dat=ClassPatternData$surfaces) {

  #--------------------------------------------------------------
  # 
  # TITLE:     surfplot()
  # AUTHOR:    TARMO REMMEL
  # DATE:      23 JULY 2013
  # CALLS:     NA
  # CALLED BY: NA
  # NEEDS:     MATRIX OBJECT FOR SURFACE AND X,Y COORDINATES
  #            TO DEFINE THE PROPORTION AND RHO VALUES TO PLOT
  #            REQUIRES surfaces OBJECT AS LOOKUP TABLE
  #            REQUIRES AN INTEGER FOR metric TO INDICATE WITH ONE TO WORK WITH
  # NOTES:     prop must be > 0.12
  #            rho must be > 0.1 
  #            IF colour=FALSE, THE PLOT IS DONE IN BW
  #            IF drop=FALSE, THE DROP LINE FROM THE POINT IS OMITTED
  #--------------------------------------------------------------

  surfaces <- dat

  plot.new()
  if(cross) {
    par(pty="s", mfrow=c(1,3))
  }
  else {
    par(pty="s", mfrow=c(1,1))
  }
  
  # PLOT PERSPECTIVE SURFACE WITH PROPORTION AND RHO POINT INDICATED WITH DROP LINE
  if(colour) {
    surfaceobj <- apply(surfaces[metric,,,], MARGIN=c(1,2), median)
    surf <- persp(seq(0.1,0.9,by=0.1), seq(0,0.2499999, by=0.2499999/10)*4, surfaceobj, ticktype="detailed", xlab="Proportion", ylab="Rho", zlab="Metric", theta=-45)
    if(drop) {
      from <- trans3d(x=prop, y=rho, z=surfaceobj[(prop*9)+1, (rho*11)+1], surf)
      to <- trans3d(x=prop, y=rho, z=min(surfaceobj), surf)
      segments(from$x, from$y, to$x, to$y, col="Red", lwd=1, lty="dotted")
    }
    points(trans3d(x=prop, y=rho, z=surfaceobj[(prop*9)+1, (rho*11)+1], surf), col="Red", pch=19)
  }
  else {
    surfaceobj <- apply(surfaces[metric,,,], MARGIN=c(1,2), median)
    surf <- persp(seq(0.1,0.9,by=0.1), seq(0,0.2499999, by=0.2499999/10)*4, surfaceobj, ticktype="detailed", xlab="Proportion", ylab="Rho", zlab="Metric", theta=-45)
    if(drop) {
      from <- trans3d(x=prop, y=rho, z=surfaceobj[(prop*9)+1, (rho*11)+1], surf)
      to <- trans3d(x=prop, y=rho, z=min(surfaceobj), surf)
      segments(from$x, from$y, to$x, to$y, lwd=1, lty="dotted")
    }
    points(trans3d(x=prop, y=rho, z=surfaceobj[(prop*9)+1, (rho*11)+1], surf), pch=19)
  
  }
  
  if(cross) {
    # PLOT BOXPLOTS ACROSS THE 11 LEVELS OF SPATIAL AUTOCORRELATION (RHO)
    plot(factor(round(seq(0,0.2499999, by=0.2499999/10)*4, 2)), surfaces[metric,round(rho*11),,], xlab="Spatial Autocorrelation", ylab="Metric Value")
    title("Metric Versus Autocorrelation (Rho)")

    # PLOT BOXPLOTS ACROSS THE 9 LEVELS OF PROPORTION
    plot(factor(seq(0.1,0.9,by=0.1)), surfaces[metric,,round(prop*9),], xlab="Proportion", ylab="Metric Value")
    title("Metric Versus Proportion")
  }
    
}
