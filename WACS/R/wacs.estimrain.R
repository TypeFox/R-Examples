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

wacs.estimrain=function(DATA,rain.model="Gamma",method="MLE",plot.it=FALSE,DIR="./"){
  #####################################################
  #
  #
  # wacs.estimrain: estimates the parameters the rain distribution, 
  #                 creates normal scores and does some control plots
  #
  # ARGUMENTS 
  #     DATA        : array of DATA, as produced by the function WACSdata
  #     rain.model  : "Gamma" for a Gamma distribution; 
  #                      "AB" for the Allard-Bourotte model 
  #                    "None" if no transformation
  #     method      : "MLE" for Maximul Likelihood 
  #                   "MOM" for Method of Moment
  #     plot.it     : whether plots are produced or not
  #     DIR         : directory in which plots are printed
  #
  # VALUE
  #    a list of parameters and the normal scores
  #
  # CALLS
  #    estim.Gamma
  #####################################################
  
  Ns   = max(DATA$season)
  Nd   = dim(DATA)[1]
  
  if ( (rain.model != "Gamma") & (rain.model != "None")) stop("[wacs.estimrain] rain.model should 'Gamma' or 'None'")
  # Estimating the parameters of the Gamma transform 
  
  if(rain.model=="Gamma"){ 
    Nparam   = 2
    par.rain = array(-999,c(Ns,Nparam))
    for (s in 1:Ns){
      par.rain[s,] = estim.Gamma(DATA$rain[(DATA$season == s) & (DATA$rain > 0)],method=method,seas=s,plot.it=plot.it,DIR=DIR)
    }  
    # Creating the normal score transform of rain  
    NscTR=rep(-999,Nd)
    for (j in 1:Nd){
      if (DATA$rain[j] > 0){
        s = DATA$season[j]
        NscTR[j] = qnorm(pgamma(DATA$rain[j],scale=par.rain[s,1],shape=par.rain[s,2]))
      }
    }
  }else{}
  
  if(rain.model=="None"){
    Nparam = 0
    par.rain = NULL
    NscTR=rep(-999,Nd)
    for (j in 1:Nd){
      if (DATA$rain[j] > 0){
        s = DATA$season[j]
        NscTR[j] = DATA$rain[j]
      }
    }
  }else{}
  
  # Estimating the parameters of the Allard-Bourotte model 
#  if(rain.model=="AB"){ 
#    Nparam = 4
#    par.rain = array(-999,Ns,Nparam)
#    colnames(par.rain) = c("threshold","scale","shape","power")
#    par.rain[s,] = estim.AB(DATA$rain[(DATA$season ==s) & (DATA$rain > 0)],method=method,seas=s,plot.it=plot.it,DIR=DIR)
#  }else{}
  
  return(list(par.rain,NscTR))
  
  
}





estim.Gamma=function(y,method="MLE",seas,plot.it=FALSE,DIR="./"){
  #####################################################
  #
  #
  # estima.Gamma: estimates the parameters of a Gamma distribution and does 
  #                        some plots
  #
  # ARGUMENTS 
  #                 y    : vector of values
  #                 seas : season (usefull for plot.it.)
  #
  #
  # VALUE
  #    a list containing scale and shape parameters
  #####################################################
  
  if( (method != "MLE" ) & (method != "MOM") ) stop("[AdjustGamma] method must be MLE or MOM")
  if (is.logical(plot.it) == FALSE) stop("[AdjustGamma] plot.it must be TRUE or FALSE")
  
  # MLE of scale and shape parameter 
  if(method=="MLE"){
    m  = mean(y,na.rm=TRUE)
    lm = mean(log(y),na.rm=TRUE)
    
    phi = function(x,mean,logmean){
      z = abs(log(x) - digamma(x) + lm - log(m))
      return(z)
    }
    shape = optimize(phi,interval=c(0,100),mean=m,logmean=lm)$minimum
    scale = m/shape
  }else{}
  
  # MOM of scale and shape parameter 
  if(method=="MOM"){
    m = mean(y,na.rm=TRUE)
    v = var(y,na.rm=TRUE)
    
    scale = v/m
    shape = m/scale
  }else{}
  
  # Plots
  
  if (plot.it==TRUE){
    pdf(paste(DIR,"Precipitation_season",seas,".pdf",sep=""))
    par(mfrow=c(1,2))
    
    plot(sort(y[!is.na(y)]),qgamma(c(1:length(y[!is.na(y)]))/length(y[!is.na(y)]),
                                   scale=scale,shape=shape),type="l",ylab="Observed quantiles",
         xlab="Theoretical quantiles",main=paste("Gamma qq-plot season=",seas,sep=""))
    abline(0,1,col="blue")
    legend("topleft",legend=c(paste("scale=",round(scale,2)),paste("shape=",round(shape,2))))
    
    hist(y[y!=0],freq=FALSE,xlab="Precipitation intensity",ylab="density",main=paste("season=",seas,sep=""),
         n=round(max(y[!is.na(y)])+1,0))
    p=seq(0,max(y[!is.na(y)]),by=0.1)
    d=dgamma(p,scale=scale,shape=shape)
    lines(p,d,col="blue",lwd=3)
    dev.off()
  }else{}
  return(c(scale,shape))
}





