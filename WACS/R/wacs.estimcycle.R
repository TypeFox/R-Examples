  ###################################################################
  #
  # These functions are part of WACSgen V1.0 
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
  
  wacs.estimcycle = function(y,spar=0.7,trend.norm="L2",plot.it=TRUE,DIR="./"){
  #####################################################
  #
  # wacs.estimcycle  : removes the seasonal component of the values. 
  #                    returns the residuals, the central tendency 
  #                    and a measure of deviation
  #
  # ARGUMENT 
  #     y         :    vector of values; length must be a multiple of 365
  #     spar      :    smoothness paramater for seasonal trend
  #     trend.norm:    norm to be used 
  #                    "L2" = quadratic norm, leading to mean and st
  #                    "L1" = Absolute value norm, leading to median and absolute deviation
  #     plot.it   : whether plots are produced or not
  #     DIR       : directory in which plots are printed
  #
  # VALUE
  #     A list containing three elements: 
  #     1/ an array containing the residuals (Nrow = Nd, number of days; 
  #                                           Ncol = Nv, number of variables to be centered) 
  #     2/ 365 x Nv array of central tendencies. 
  #     3/ 365 x Nv array of central deviations. 
  #
  # CALLS
  #     extract.annual.trend
  #####################################################

  Central   = NULL
  Deviation = NULL
  Resid     = NULL
  Nv        = dim(y)[2]
  for (v in 1:Nv)
  {
    M     = matrix(y[,v],ncol=365,byrow=TRUE)
    Dv    = extract.annual.trend(y[,v],spar=spar,trend.norm=trend.norm)
    Central   = cbind(Central,Dv[,1])
    Deviation = cbind(Deviation,Dv[,2])
    Resid     = cbind(Resid,(y[,v] - Dv[,1])/Dv[,2])
    colnames(Resid)[v] = paste("V",v,"d",sep="")
    
    if (plot.it==TRUE){
      png(paste(DIR,"Annual_",colnames(y)[v],".png",sep=""))
      boxplot(t(M)~c(1:365),outline=FALSE,xlab="Julian day",ylab=paste("Variable V",v,sep=""),na.action=drop)
      lines(Dv[,1],col="red",lwd=3)
      lines(Dv[,1]-2*Dv[,2],col="red",lwd=1)
      lines(Dv[,1]+2*Dv[,2],col="red",lwd=1)
      dev.off()
    }
  }
  return(list(Resid,Central,Deviation))
}


extract.annual.trend = function(y,spar=0.7,trend.norm="L2"){
  #####################################################
  #
  # extract.annual.trend :extract the central tendency and a measure of deviation
  #
  # ARGUMENTS 
  #                  y   : vector of values; must be of length 365
  #                  spar: smoothness paramater for seasonal trend
  #            trend.norm: norm to be used 
  #                        "L2" = quadratic norm, leading to mean and st
  #                        "L1" = Absolute value norm, leading to median and absolute deviation
  #
  # VALUE
  #     a 365 x 2 array. First column:  vector of central values; 
  #                      Second column: vector of deviation    
  #
  #####################################################
  M=matrix(y,ncol=365,byrow=TRUE)
  if (trend.norm=="L1"){
    y.center = apply(M,2,median,na.rm=TRUE)
    dev = function(v) z=mean(abs(v-median(v,na.rm=TRUE)),na.rm=TRUE)
  }else{}
  if (trend.norm=="L2"){
    y.center = apply(M,2,mean,na.rm=TRUE)
    dev = function(v) z = sqrt(mean( (v-mean(v,na.rm=TRUE))^2, na.rm=TRUE))
  }else{}
  y.dev = apply(M,2,dev)
  
  # smoothing
  y.center.smooth = smooth.spline(rep(y.center,3),spar=spar)[[2]][366:730]
  y.dev.smooth    = smooth.spline(rep(y.dev,3),spar=spar)[[2]][366:730]
  return(cbind(y.center.smooth,y.dev.smooth))
}
  
  
  
