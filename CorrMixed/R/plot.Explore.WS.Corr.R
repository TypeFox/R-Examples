plot.Explore.WS.Corr <- function(x, Est.Corrs=TRUE, Indiv.Corrs=FALSE, Add.CI=FALSE, 
                                 Add.CI.Smoothed=TRUE, Smoother.Span=0.2, 
                                 Add.Boot.Corrs=FALSE, Add.CI.Polygon=FALSE, ylim=c(-1, 1), 
                                 xlab="Time Lag", ylab="Reliability", ...){
  
  F_val <- Smoother.Span
  
  if (missing(xlab)) {xlab="Time lag"}
  if (missing(ylab)) {ylab="Reliability"}   
  
  Object <- x 
  
  plot(x=1:length(Object$Est.Corr), y=Object$Est.Corr, col=0, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  if (Indiv.Corrs==TRUE) {points(x = Object$All.Corrs[,1], 
                                 y=Object$All.Corrs[,2], col="grey", pch=19, cex=1)}
  
  
  if (Add.Boot.Corrs==TRUE){
    for (i in 1: dim(Object$Bootstrapped.Corrs)[1]){
      lines(x = 1:length(Object$Bootstrapped.Corrs[i,]),
            y=Object$Bootstrapped.Corrs[i,], col="grey")
    }
  }
  
  if (Est.Corrs==TRUE){
    lines(y=Object$Est.Corr, x=1:length(Object$Est.Corr), lwd=2, col=1)
  }
  
  if (Add.CI==TRUE){
    lines(x= (1: length(Object$CI.Lower)), y=Object$CI.Lower, lty=2, lwd=2)  
    lines(x= (1: length(Object$CI.Upper)), y=Object$CI.Upper, lty=2, lwd=2)  
  }
  
  
  if (Add.CI.Smoothed==TRUE){
    plot_low <- lowess(x = (1: length(Object$CI.Lower)), y=Object$CI.Lower, f=F_val)$y 
    plot_high <- lowess(x = (1: length(Object$CI.Upper)), y=Object$CI.Upper, f=F_val)$y
    lines(y=plot_low, x=(1:length(plot_low)), lty=2, lwd=2)  
    lines(y=plot_high, x=(1:length(plot_high)), lty=2, lwd=2)    
  }
  
  
  if (Add.CI.Polygon==TRUE){
    plot_low <- lowess(x = (1: length(Object$CI.Lower)), y=Object$CI.Lower, f=F_val)$y 
    plot_high <- lowess(x = (1: length(Object$CI.Upper)), y=Object$CI.Upper, f=F_val)$y
    lines(y= plot_low, x=(1:length(plot_low)), col="ivory2")
    lines(y= plot_high, x=(1:length(plot_high)), col="ivory2")
    
    polygon(c(1:length(plot_low), rev(1:length(plot_low))), c(plot_high, rev(plot_low)), 
            col = "ivory2", border = NA)
  }
}
