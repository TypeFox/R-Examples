###################################################################
#
# This function is part of WACSgen V1.0 
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trépos MIA-T, INRA
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################

#' Produces validation and/or WACS comparison plots
#' 
#' For plotting validation figures from outputs generated when calling \link{WACSvalid} or \link{WACScompare}.
#' Figures are either displayed or printed in a file
#'
#' @export
#' 
#' @param    wacsvalid   Output, as obtained when calling \link{WACSvalid} or \link{WACScompare}
#' @param    file        File in which to write the figure. Default is \code{NULL}. If \code{file=NULL}, no file 
#'                       is produced; the figure is only produced in the graphical interface 
#' 
#' @return  No return. A Figure is either displayed or printed in a file.
#' 
#' @examples
#' \dontrun{
#'   ## Simple example
#'   data(ClimateSeries)
#'   ThisData = WACSdata(ClimateSeries,from="1995-01-01",to="1999-12-31")
#'   ThisPar  = WACSestim(ThisData)
#'   ThisSim  = WACSsimul(ThisPar,from="1995-01-01",to="1999-12-31")
#'   ThisVal  = WACSvalid(what="Sim",wacsdata = ThisData, wacspar = ThisPar, 
#'                        wacssimul = ThisSim,varname="tmin")
#'   WACSplot(ThisVal,file="ThisFile.pdf")
#'}
#' 

WACSplot = function(wacsvalid,
                    file=NULL)
{
  #  Some checking
  
  if (!(class(wacsvalid) %in% c("WACSvalidSim","WACSvalidRain","WACSvalidMeanSd","WACSvalidBiVar",
                                "WACSvalidCorTemp","WACSvalidSumBase","WACSvalidPersistence"))){
    stop("[WACSplot] wacsvalid should be of a correct class")
  }
  
  Lab1 = wacsvalid$labels[1]
  Lab2 = wacsvalid$labels[2]
  #
  # Different types of classes
  #
  
  if (class(wacsvalid) == "WACSvalidSim") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(1,1))
    vmax = max(wacsvalid$zObs,wacsvalid$zSim,wacsvalid$CentralObs+2*wacsvalid$DeviationObs,
               wacsvalid$CentralSim+2*wacsvalid$DeviationSim)
    vmin = min(wacsvalid$zObs,wacsvalid$zSim,wacsvalid$CentralObs-2*wacsvalid$DeviationObs,
               wacsvalid$CentralSim-2*wacsvalid$DeviationSim)
    plotmax = vmax + 0.01*(vmax-vmin)
    plotmin = vmin - 0.01*(vmax-vmin)
    plot(1:365,wacsvalid$CentralObs,type="l",lwd=3,col="brown",ylim=c(plotmin,plotmax),
         xlab="Julian days",ylab=wacsvalid$varname,main=paste(Lab1, "(brown) and", Lab2, "(blue)",wacsvalid$varname))
    for (i in 1: min(dim(wacsvalid$zSim)[1], dim(wacsvalid$zObs)[1])){
      points(1:365,wacsvalid$zSim[i,],type="l",lty=1,lwd=1,col="lightblue")
      points(1:365,wacsvalid$zObs[i,],type="l",lty=1,lwd=1,col="lightsalmon")
    }
    points(1:365,wacsvalid$CentralObs,type="l",lwd=3,col="brown")
    points(1:365,wacsvalid$CentralObs+2*wacsvalid$DeviationObs,type="l",lwd=1,col="brown")
    points(1:365,wacsvalid$CentralObs-2*wacsvalid$DeviationObs,type="l",lwd=1,col="brown")
    points(1:365,wacsvalid$CentralSim,type="l",lty=2,lwd=3,col="blue")
    points(1:365,wacsvalid$CentralSim-2*wacsvalid$DeviationSim,type="l",lty=2,lwd=1,col="blue")
    points(1:365,wacsvalid$CentralSim+2*wacsvalid$DeviationSim,type="l",lty=2,lwd=1,col="blue")  
    seasons = wacsvalid$seasons
    ns = dim(seasons)[1]
    for (i in 1:ns){
      day = paste("01/",seasons[i,1],"/",seasons[i,2],sep="")
      abline(v=wacs.dayOfYear(day),lty=3,col="gray")
    }
    if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidRain") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(length(wacsvalid),2))
    for (s in 1:(length(wacsvalid)-1)) {#TODO have fixed length(wacsvalid)
      scale=wacsvalid[[s]]$par[1];
      shape=wacsvalid[[s]]$par[2];
      #qq plot
      plot(wacsvalid[[s]]$theoretical,wacsvalid[[s]]$observed,type="l",
              ylab=paste(Lab1,"quantiles"), xlab=paste(Lab2,"quantiles"),
              main=paste("Gamma qq-plot season=",s,sep=""))
      abline(0,1,col="blue")
      legend("topleft",legend=c(paste("scale=",round(scale,2)),
                      paste("shape=",round(shape,2))))
      #hist
      hist(wacsvalid[[s]]$observed,freq=FALSE,xlab="Precipitation intensity",
              ylab="density",main=paste("season=",s,sep=""),
           n=round(max(wacsvalid[[s]]$observed)+1,0))
      p=seq(0,max(wacsvalid[[s]]$observed),by=0.1)
      d=dgamma(p,scale=scale,shape=shape)
      lines(p,d,col="blue",lwd=3) 
    }
    if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidMeanSd") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(2,2));
      boxplot(wacsvalid$meanObs,col="brown",las=1,range=0,
              main=paste(Lab1, "mean,",wacsvalid$varname),
              ylim=c(min(wacsvalid$meanObs,wacsvalid$meanSim),
                      max(wacsvalid$meanObs,wacsvalid$meanSim)));
      g = boxplot(wacsvalid$meanObs,plot=F)
      
      boxplot(wacsvalid$meanSim,col="cyan",las=1,range=0,
              main=paste(Lab2, "mean,",wacsvalid$varname),
              ylim=c(min(wacsvalid$meanObs,wacsvalid$meanSim),
                      max(wacsvalid$meanObs,wacsvalid$meanSim)));
      lines(1:12,g$stat[3,],col="brown",type="l",lwd=3)
      lines(1:12,g$stat[2,],col="brown",type="l",lty=2,lwd=2)
      lines(1:12,g$stat[4,],col="brown",type="l",lty=2,lwd=2)
      
      boxplot(wacsvalid$sdObs,col="brown",las=1,range=0,
              main=paste(Lab1,"sd,",wacsvalid$varname),
              ylim=c(min(wacsvalid$sdObs,wacsvalid$sdSim),
                      max(wacsvalid$sdObs,wacsvalid$sdSim)));
      g = boxplot(wacsvalid$sdObs,plot=F)
      
      boxplot(wacsvalid$sdSim,col="cyan",las=1,range=0,
              main=paste(Lab2,"sd,",wacsvalid$varname),
              ylim=c(min(wacsvalid$sdObs,wacsvalid$sdSim),
                      max(wacsvalid$sdObs,wacsvalid$sdSim)));
      lines(1:12,g$stat[3,],col="brown",type="l",lwd=3)
    lines(1:12,g$stat[2,],col="brown",type="l",lty=2,lwd=2)
    lines(1:12,g$stat[4,],col="brown",type="l",lty=2,lwd=2)
    if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidCorTemp") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(1,2));
      boxplot(wacsvalid$corObs,
              main=paste(Lab1,"Temp cor,",wacsvalid$varname),
              ylim=c(min(wacsvalid$corObs,wacsvalid$corSim,na.rm=TRUE),
                      max(wacsvalid$corObs,wacsvalid$corSim,na.rm=TRUE)),
              range=0,col="brown");
      h = boxplot(wacsvalid$corObs,plot=F);
      boxplot(wacsvalid$corSim,
              main=paste(Lab2,"Temp Cor,",wacsvalid$varname),
              ylim=c(min(wacsvalid$corObs,wacsvalid$corSim,na.rm=TRUE),
                      max(wacsvalid$corObs,wacsvalid$corSim,na.rm=TRUE))
              ,range=0,col="cyan")
      lines(1:12,h$stat[3,],col="brown",type="l",lwd=3)
      lines(1:12,h$stat[2,],col="brown",type="l",lty=2,lwd=2)
      lines(1:12,h$stat[4,],col="brown",type="l",lty=2,lwd=2)
      if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidBiVar") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(1,2));
    boxplot(wacsvalid$corObs,
            main=paste(Lab1,"Correl,",wacsvalid$varname,",",wacsvalid$varname2),
            ylim=c(min(wacsvalid$corObs,wacsvalid$corSim),
                   max(wacsvalid$corObs,wacsvalid$corSim)),
            range=0,col="brown");
    h = boxplot(wacsvalid$corObs,plot=F);
    boxplot(wacsvalid$corSim,
            main=paste(Lab2, "Correl,",wacsvalid$varname,",", wacsvalid$varname2),
            ylim=c(min(wacsvalid$corObs,wacsvalid$corSim),
                   max(wacsvalid$corObs,wacsvalid$corSim))
            ,range=0,col="cyan")
    lines(1:12,h$stat[3,],col="brown",type="l",lwd=3)
    lines(1:12,h$stat[2,],col="brown",type="l",lty=2,lwd=2)
    lines(1:12,h$stat[4,],col="brown",type="l",lty=2,lwd=2)
    if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidSumBase") {
    if(!is.null(file)){pdf(file)}
    boxplot(cbind(wacsvalid$SumObs,wacsvalid$SumSim),
            main=paste("Sum of",wacsvalid$varname,", base",wacsvalid$base),
            ylim = c(min(wacsvalid$SumObs,wacsvalid$SumSim), 
                     max(wacsvalid$SumObs,wacsvalid$SumSim)),
            range=0,col=c("brown","blue"),names=c(Lab1,Lab2))
    h = boxplot(wacsvalid$SumObs,plot=F)$stats[3,1]
    abline(h=h,col="brown")
    if(!is.null(file)){dev.off()}
  }
  
  if (class(wacsvalid) == "WACSvalidPersistence") {
    if(!is.null(file)){pdf(file)}
    par(mfrow=c(1,1))
    ll = length(wacsvalid$FreqObs)
    Nobs = sum(wacsvalid$FreqObs*(1:ll))
    Nsim = sum(wacsvalid$FreqSim*(1:ll))
    if (wacsvalid$above){
      ch.above="above"
    }else{
      ch.above="below"
    }
    barplot(cbind(wacsvalid$FreqObs,wacsvalid$FreqSim),names.arg=c(1:ll,1:ll),beside=TRUE,col=c(rep("brown",ll),rep("blue",ll)),
            main = paste("Persistence of ",wacsvalid$varname,ch.above," base",wacsvalid$base))
    legend("topright",legend=c(paste(Lab1, "N=",Nobs),paste(Lab2,"N=",Nsim)),fill=c("brown","blue"))
    if(!is.null(file)){dev.off()}
  }
  
  
  # Closing WACSplot  
}
