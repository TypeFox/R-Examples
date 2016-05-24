`DistributionPlotGamma` <-
function(shape=1,rate=1,length=100, xlab="x", ylab="Density",
           main=bquote(paste("Gamma Distribution: ",alpha,"=",
             .(signif(shape,digits=signif.digits)),",",beta,"=",
             .(signif(rate,digits=signif.digits)))), signif.digits=3)
{
  ## Purpose: plot gamma distribution (guts borrowed from Rcmdr)
  ##
  ## shape = alpha
  ## rate = 1/scale = beta
  ## length:   no. of points to be used for plotting smooth distribution
  ## main:     title for plot
  ## ylab:     y-axis label
  ## xlab:     x-axis label
  ## signif:   number of significant digits for default 'main' title

  min <- round(qgamma(.0005, shape=shape, rate=rate))
  max <- round(qgamma(.9995, shape=shape, rate=rate))
  .x <- seq(min, max, length=length)
  plot(.x, dgamma(.x, shape=shape, rate=rate), xlab=xlab, ylab=ylab,
       main=main, type="l")
  abline(h=0, col="gray")
}

