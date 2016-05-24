`DistributionPlotNorm` <-
function(mean=0, sd=1, length=100, xlab="x", ylab="Density",
           main=bquote(paste("Normal Distribution: ",mu,"=",
             .(signif(mean,digits=signif.digits)),",",sigma,"=",
             .(signif(sd,digits=signif.digits)))), signif.digits=3)
{
  ## Purpose: plot normal distribution (guts borrowed from Rcmdr)
  ##
  ## mean:     mean of normal distribution
  ## sd:       standard deviation of normal distribution
  ## length:   no. of points to be used for plotting smooth distribution
  ## main:     title for plot
  ## ylab:     y-axis label
  ## xlab:     x-axis label
  ## signif:   number of significant digits for default 'main' title

  min <- round(qnorm(.0005, mean=mean, sd=sd), 3)
  max <- round(qnorm(.9995, mean=mean, sd=sd), 3)
  .x <- seq(min, max, length=length)
  plot(.x, dnorm(.x, mean=mean, sd=sd), xlab=xlab, ylab=ylab,
       main=main,type="l")
  abline(h=0, col="gray")
}

