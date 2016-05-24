# Name   : plot.R
# Desc   : A tweaked "plot" function designed to easily plot R objects from
#          any of the supported estimation methods.
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

plotfit.R0.R <- function#Plot the fit of a model to epidemic data
### Plots the fit of a model to epidemic data
##details<< For internal use. Called by plotfit.
##keyword<< internal

(x, ##<< Result of est.R (class R0.R)
 all=TRUE, ##<< Should the whole epidemic curve be shown
 xscale="w", ##<< Scale to be adjusted on X axis. Can be "d" (day), "w" (week (default)), "f" (fornight), "m" (month).
 SB.dist=TRUE, ##<< Should R distribution throughout the epidemic be plotted for SB method? (default: TRUE)
 ... ##<< Parameters passed to plot
 ) 
  ##details<< For internal use. Called by plotfit.R0.sR.
  
  
  # Code
  
{
  #Make sure x is of class "R0.R"
  if (class(x)!="R0.R") stop("'x' must be of class R0.R")
  
  if (x$method.code %in% c("EG","ML","TD")) {
    do.call(plotfitRxx, args=list(x=x, xscale=xscale, ...) )
  }
  else {
    do.call(paste("plotfitR",x$method.code,sep=""), args=list(x=x, xscale=xscale, SB.dist=SB.dist, ...) )
  }
  ### Called for its side effect :
  ### Draws the fit of one estimation method to the data.
}


plotfitRxx <- function(x, xscale,...)#Internal plotfit method for EG, ML and TD estimates
###Internal plotfit method for EG, ML and TD estimates
{
  ##keyword<< internal
  epid = x$epid
  
  #Get data used for the fit
  begin = x$begin
  begin.nb = x$begin.nb
  end = x$end
  end.nb = x$end.nb
  
  epid.used.for.fit = list(incid=epid$incid[begin.nb:end.nb], t=epid$t[begin.nb:end.nb])
  
  #Plot the whole original epidemic data
  plot(epid$t,epid$incid, xlab="Time",ylab="Incidence",t='s', xaxt="n", main=paste("Epidemic curve & model (", x$method,")"))
  
  #Add a line showing predicted simulation
  lines(epid.used.for.fit$t,x$pred,col='red')
  
  #Highlight the original points
  points(epid.used.for.fit$t,epid.used.for.fit$incid,pch=19)

  #Finally, X-Axis labels
  div = get.scale(xscale)
  #Where should labels be on axis
  atLab = pretty(epid$t, n=length(epid$t)/div)
  #What should labels say
  lab = format(pretty(epid$t, n=length(epid$t)/div))
  axis(1, at=atLab, labels=lab)
  
}

plotfitRAR <- function(x, xscale,...)#Internal plotfit method for AR estimates
###Internal plotfit method for AR estimates
{
  ##keyword<< internal
  epid <- x$epid
  epid.orig <- x$epid.orig
  epid.used.for.fit = list(incid=epid.orig$incid, t=epid.orig$t)
  
  #Plot the whole original epidemic data
  plot(epid$t,epid$incid, xlab="Time",ylab="Incidence",t='s', xaxt="n", main="Epidemic curve (Attack Rate)")
  
  #Highlight the original points
  points(epid.used.for.fit$t,epid.used.for.fit$incid,pch=19)
  
  #Finally, X-Axis labels
  div = get.scale(xscale)
  #Where should labels be on axis
  atLab = pretty(epid$t, n=length(epid$t)/div)
  #What should labels say
  lab = format(pretty(epid$t, n=length(epid$t)/div))
  axis(1, at=atLab, labels=lab)
}

plotfitRSB <- function(x, xscale, SB.dist,...)#Internal plot method for SB estimates
###Internal plotfit method for SB estimates
{
  ##keyword<< internal
  epid = x$epid
  
  #Get data used for the fit
  begin = x$begin
  begin.nb = x$begin.nb
  end = x$end
  end.nb = x$end.nb
  
  epid.used.for.fit = list(incid=epid$incid[begin.nb:end.nb], t=epid$t[begin.nb:end.nb])
  
  #Plot the whole original epidemic data
  plot(epid$t,epid$incid, xlab="Time",ylab="Incidence",t='s', xaxt="n", main=paste("Epidemic curve & model (", x$method,")"))
  
  #Add a line showing predicted simulation
  lines(epid.used.for.fit$t,x$pred,col='red')
  
  #Highlight the original points
  points(epid.used.for.fit$t,epid.used.for.fit$incid,pch=19)
  
  #Finally, X-Axis labels
  div = get.scale(xscale)
  #Where should labels be on axis
  atLab = pretty(epid$t, n=length(epid$t)/div)
  #What should labels say
  lab = format(pretty(epid$t, n=length(epid$t)/div))
  axis(1, at=atLab, labels=lab)
  
  #When plotting Bayesian, if SB.dist is enabled, plot some R distributions throughout the epidemic
  if (SB.dist == TRUE) {
    #x11()
    dev.new()
    split.screen(c(3,3))
    if (end.nb-begin.nb>8) {
      num.to.plot <- c(1, rep(NA, 8))
    }
    else {
      num.to.plot <- c(begin.nb:end.nb)
    }
    for (i in 1:length(num.to.plot)) {
      if (i == 1) {
        screen(1)
        plot(y=x$proba.Rt[[num.to.plot[i]]], x=seq(from=0, to=(length(x$proba.Rt[[num.to.plot[i]]])/100-0.01), by=0.01), xlab="R value", ylab="PDF", type="l", main=paste("t=",num.to.plot[i]))
        abline(v=(which.max((cumsum(x$proba.Rt[[num.to.plot[i]]])) >= 0.025)-1)/100, col="red", lty="dotted")
        abline(v=(which.max((cumsum(x$proba.Rt[[num.to.plot[i]]])) >= 0.975)-1)/100, col="red", lty="dotted")
        next
      }
      if (is.na(num.to.plot[i])) {
        num.to.plot[i] = num.to.plot[i-1] + floor(length(x$epid$incid[begin.nb:end.nb])/9)
      }
      
      screen(i)
      plot(x$proba.Rt[[num.to.plot[i]]], x=seq(from=0, to=(length(x$proba.Rt[[num.to.plot[i]]])/100-0.01), by=0.01), xlim=c(0,((length(x$proba.Rt[[num.to.plot[i]]]) - which.max(rev(x$proba.Rt[[num.to.plot[i]]])>0) + 1))/100 - 0.01), xlab="R value", ylab="PDF", pch=NA_integer_, type="l", main=paste("t=",num.to.plot[i]))
      abline(v=(which.max((cumsum(x$proba.Rt[[num.to.plot[i]]])) >= 0.025)-1)/100, col="red", lty="dotted")
      abline(v=(which.max((cumsum(x$proba.Rt[[num.to.plot[i]]])) >= 0.975)-1)/100, col="red", lty="dotted")
      
    }
    #Closing devices
    close.screen(all.screens=TRUE)
  }
}
