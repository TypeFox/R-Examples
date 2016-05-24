ffplot <- function (t,y,skeleton,fcast,std.error) {
  #*********************************************************************
  #  This function makes a plot of the flow field forecast
  #
  #  Input: t - time series observation times
  #         y - time series response values
  #         skeleton - Matrix containing the data skeleton
  #         fcast - vector holding the forecast values
  #         std.error - vector holding the standard errors
  #
  #  Output: Plot of times series and forecast values
  #
  #  References: None   
  #
  #**********************************************************************
  par(mai=c(1.1,1.0,0.9,1.0))
  
  knots <- length(skeleton$sd)-4
  space <- skeleton[knots+2,2]
  lknot <- skeleton[knots+3,2]
  fknot <- skeleton[knots+3,1]
  steps <- length(fcast)
  
  kvector <- seq(fknot-space,lknot,space)
  tf <- rep(0,steps) 
  for (i in 1:steps){
    tf[i] <-lknot + i*steps 
  }
  sd <- skeleton$sd[1:(knots+1)]
  xmin <- (fknot - space)
  xmax <- lknot + steps*space
  ymin <- min(y)  
  ymax <- max(y) 

  # Plot the raw time series data
  plot(t,y,pch=20,bg="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='L',main="Flow Field Forecast")

  # Plot the forecast values
  points(tf,fcast,pch=21,col="blue")

  # Plot the penalized spline regression
  points(kvector,sd,type="l",col="red")
  
  # Plot the error bands
  points(tf,fcast+2*std.error,type="l",col="green")
  points(tf,fcast-2*std.error,type="l",col="green")
  
  legend("topright", inset=c(0,0),c("TS Data","Forecast","PSR","+/- 2SE"),cex=0.8,col=c("black","blue","red","green"),bty="n",pch=c(16,1,NA,NA),lty=c(NA,NA,1,1))
}