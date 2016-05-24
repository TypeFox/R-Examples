# Name   : plot.R
# Desc   : A tweaked "plot" function designed to easily plot R objects from
#          any of the supported estimation methods.
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

plot.R0.R <- function#Plot the R0/Rt value along with confidence interval
###Plot the R0/Rt value along with confidence interval
##details<< For internal use. Called by plot.
##keyword<< internal

(x, ##<< Result of est.R (class R)
 all=TRUE, ##<< Should the whole epidemic curve be shown
 xscale="w", ##<< Scale to be adjusted on X axis. Can be "d" (day), "w" (week (default)), "f" (fornight), "m" (month).
 TD.split=FALSE, ##<< Parameter to force the display of both R(t) and the epidemic curve in the same window for TD method
 ... ##<< Parameters passed to plot
) 

  
# Code
  
{
	#Make sure x is of class "R0.R"
	if (class(x)!="R0.R") stop("'x' must be of class R0.R")
  
  do.call(paste("plotR",x$method.code,sep=""), args=list(x=x, xscale=xscale, TD.split=TD.split, ...) )

### Called for its side effect :
### Draws all R0 or R(t) value from one estimation methods.
}


# plot for EG : 
# this plot the exponential growth parameter r, along with R
plotREG <- function(x, ...)#Internal plot method for EG estimates
###Internal plot method for EG estimates
  
{
  ##keyword<< internal
    R <- x$R
    conf.int <- x$conf.int
    r <- x$r
    
    r.seq = seq(from=-1, to=1, by=0.01)
    exp.r.seq = exp(r.seq*x$GT$mean)
    constant.r.seq = 1+r.seq*x$GT$mean
    
    plot(y=exp.r.seq, x=r.seq, 
         xlim=c(-1,1), ylim=c(0,3), 
         type="l", lty="dashed", 
         xlab="r (Growth rate)", ylab="R (Reproduction number)", 
         main=paste("Reproduction number (", x$method,")"))
    lines(y=constant.r.seq, x=r.seq, xlim=c(-1,1), ylim=c(0,3), lty="solid")
    
    points(y=R, x=r, col="blue")
    arrows(r,conf.int[2], r, conf.int[1], angle=90, code=3, col="blue", length=0.03)
    legend(x="bottomright", c("exp(r*mean(GT))", "1+r*mean(GT)", "R"), pch=c(NA_integer_, NA_integer_, 21), 
      col=c("black", "black", "blue"), lty=c("dashed", "solid", "blank"), merge=TRUE)
}

plotRAR <- function(x, ...)#Internal plot method for AR estimates
###Internal plot method for AR estimates
{
  ##keyword<< internal
   	R <- x$R
    conf.int <- x$conf.int
    
    plot(x=1, y=R, col="blue", 
         ylab="R0 value", xaxt="n", 
         xlab="", 
         main=paste("Reproduction number (", x$method,")"))
    arrows(1,conf.int[2], 1, conf.int[1], angle=90, code=3, col="blue", length=0.03)   
}


plotRML <- function(x, ...)#Internal plot method for ML estimates
###Internal plot method for ML estimates
{
  ##keyword<< internal
   	R <- x$R
    conf.int <- x$conf.int
    
    plot(x=1, y=R, col="blue", 
         ylab="R (Reproduction number)", xaxt="n", xlab="", 
         main=paste("Reproduction number (", x$method,")"))
    arrows(1,conf.int[2], 1, conf.int[1], angle=90, code=3, col="blue", length=0.03)   
}

plotRTD <- function(x, xscale, TD.split,...)#Internal plot method for TD estimates
###Internal plot method for TD estimates
{
  ##keyword<< internal
    #We plot R(t) values along with 95CI
    R = x$R
    epid = x$epid
    epid$t = epid$t[x$begin.nb:(x$end.nb-1)]
    polygon.x = c(epid$t, rev(epid$t))
    polygon.y = c(x$conf.int[1:length(epid$t),1], rev(x$conf.int[1:length(epid$t),2]))
    
    div = get.scale(xscale)
    #Where should labels be on axis
    atLab = pretty(epid$t, n=length(epid$t)/div)
    #What should labels say
    lab = format(pretty(epid$t, n=length(epid$t)/div))

    if (TD.split == TRUE) {
      #par(bg="white")
      split.screen(c(2,1))
      screen(1)
      plot(epid$incid[x$begin.nb:x$end.nb], 
           type='s', xlab="Time", ylab="Incidence", 
           xaxt="n", main="")
      
      axis(1, at=atLab, labels=lab)
      screen(2)
    }
    plot(epid$t, R[1:length(epid$t)], ylim=c(0, max(polygon.y)), xlab="Time", ylab="R(t)", xaxt="n", pch=NA_integer_, lty="blank", main=paste("Reproduction number (", x$method,")"))
    polygon(polygon.x, polygon.y, col="gray", border=NA)
    lines(epid$t, R[1:length(epid$t)])
    abline(h=1, lty="dashed", col="gray40")
    
    #Blue = lowest quantile (default: 5%)
    #Red = highest quantile (default: 95%)
    #points(epid$t, Rt.quant$CI.lower.[1:length(epid$t)], col="blue", xaxt="n", pch=NA_integer_)
    #points(epid$t, Rt.quant$CI.upper.[1:length(epid$t)], col="red", xaxt="n", pch=NA_integer_)
    axis(1, at=atLab, labels=lab)
    
    if (TD.split == TRUE) {
      close.screen(split.screen())
    }
}

plotRSB <- function(x, xscale,...)#Internal plot method for SB estimates
###Internal plot method for SB estimates
{
  ##keyword<< internal
	  #We plot R(t) values along with 95CI
	  R= x$R
	  epid = x$epid
	  
	  epid$t = epid$t[x$begin.nb:(x$end.nb-1)]
	  polygon.x = c(epid$t, rev(epid$t))
	  polygon.y = c(x$conf.int[1:length(epid$t),1], rev(x$conf.int[1:length(epid$t),2]))
	  
	  #plot(epid$t, Rt.quant$R.t.[1:length(epid$t)], ylim=c(0, max(polygon.y)), xlab="Time", ylab="R(t)", xaxt="n", main=paste("Reproduction number (", x$method,")"),...)
	  plot(epid$t, R[1:length(epid$t)], ylim=c(0, max(polygon.y)), xlab="Time", ylab="R(t)", xaxt="n", pch=NA_integer_, lty="blank", 
    main=paste("Reproduction number (", x$method,")"))
	  polygon(polygon.x, polygon.y, col="gray", border=NA)
	  lines(epid$t, R[1:length(epid$t)])
	  abline(h=1, lty="dashed", col="gray40")

    div = get.scale(xscale)
	  #Where should labels be on axis
	  atLab = pretty(epid$t, n=length(epid$t)/div)
	  #What should labels say
	  lab = format(pretty(epid$t, n=length(epid$t)/div))
	  
    axis(1, at=atLab, labels=lab)
}


get.scale <- function(scale)#Internal scaling function to display proper X-Axis labels
###Internal scaling function to display proper X-Axis labels
{
  ##keyword<< internal
  #Scale parameters are used to adjust dates on X-Axis labels
  if (scale == "d") {
    div = 1
  } else if (scale == "w") {
    div = 7
  } else if (scale == "f") {
    div = 15
  } else if (scale == "m") {
    div = 30
  } else {
    stop("Invalid scale parameter.")
  }
  return(div)
}
