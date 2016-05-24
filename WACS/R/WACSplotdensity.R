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

#' For plotting fitted bivariate densities of residuals
#' 
#' @export
#' 
#' @param wacsdata WACS data obtained when calling \link{WACSdata} on original climate series
#' @param wacspar  WACS parameters estimated when calling \link{WACSestim}
#' @param season   season to be considered (a scalar)
#' @param dimens   a vector of length 1 or 2 indicating the marginals to be plotted
#' @param dry      indicates whether dry weather states (if \code{dry=TRUE}) or wet weather states (if \code{dry=FALSE}) must be considered
#' @param DIR      Directory in which saving the Figures
#'  
#' @note  If \code{length(dimens)=1}, the bivariate density of the variable at days \code{(d,d+1)} is plotted. 
#' If \code{length(dimens)=2}, the same-day bivariate density of the pair of variables is plotted.
#'    
#' @examples 
#' \dontrun{
#'   ## Simple example
#'   data(ClimateSeries)
#'   ThisData = WACSdata(ClimateSeries,from="1995-01-01",to="1999-12-31")
#'   ThisPar  = WACSestim(ThisData)
#'   WACSplotdensity(ThisData,ThisPar,season=2,dimens=c(2,3),dry=TRUE) 
#' }
#'

WACSplotdensity = function(wacsdata = NULL,
                           wacspar = NULL,
                           season,
                           dimens=c(1,2),
                           dry=T,
                           DIR="./"){


  
  # Initialization
  Nv = dim(wacsdata$data)[2] - 4
  Nd = dim(wacsdata$data)[1]
  Ny = Nd%/%365
  sameday = T
  if (length(dimens)==1){ 
    sameday = F}
  if (length(dimens) > 2){
    stop("WACSplotdensity: length(dimens) > 2")
  }
  if ( (1 %in% dimens) & dry){
    stop("WACSplotdensity: cannot plot rain values for DRY days")
  }
  if ( (length(dimens)==2) & (dimens[1] ==  dimens[2]) ){
    stop("WACSplotdensity: values in vector 'dimens' must be different")
  }
  s = season
  Season.par  = wacspar[[paste("Season",sep="_",s)]]
  
  # Centering and transforming rain
  Y = matrix(0,Nd,Nv)
  if(wacspar$Rain$RainModel=="Gamma"){
    par.rain = wacspar$Rain$RainPar[s,]
    Y[,1] = rep(-999,Nd)
    for (j in 1:Nd){
      if (wacsdata$data$rain[j] > 0) {
        Y[j,1] = qnorm(pgamma(wacsdata$data$rain[j],scale=par.rain[1],shape=par.rain[2]))
      }
    }
  }
  for (v in 2:Nv){ 
    Y[,v] = (wacsdata$data[,(4+v)] - rep(wacspar$Trend$Central[,v-1],Ny))/rep(wacspar$Trend$Deviation[,v-1],Ny)
  }
  
  # Stationary probability for each weather type
  TransM = Season.par$TransM
  Nt     = dim(TransM)[1]
  Prob   = abs(eigen(t(TransM))$vectors[,1])
  Prob   = Prob/sum(Prob)
  
  # List of weather types
  vdry = as.numeric(dry)
  if (dry){
    Na = 1
    Nz = Season.par$NumbWT[1]
    B  = (Y[,1] < -998) & (wacsdata$data$season == s)
  }else{
    Na = Season.par$NumbWT[1] + 1
    Nz = Season.par$NumbWT[1] + Season.par$NumbWT[2]
    B = (Y[,1]  > -998) & (wacsdata$data$season == s)
  }
  wt = Na:Nz

  # Initialization for plots
  uu = seq(-3.3,3.3,by=0.15)

  # Same variable, two succesive days 
  if (!sameday){ # One variable; day d and d+1
    V1  = dimens[1]
    varname = names(wacsdata$data)[4+V1]
    h   = wacsdensity_1d(V1,uu,wt,Season.par,Prob,vdry)
    g   = wacsdensity_time(V1,uu,wt,Season.par,Prob,vdry)
    
    Y1 = Y[B,V1] 

    if (dry){
      pdf(paste(DIR,"CSN_Season_",s,"_dry_",varname,".pdf",sep=""))
    }else{
      pdf(paste(DIR,"CSN_Season_",s,"_wet_",varname,".pdf",sep=""))
    }
    par(mfrow=c(1,2))
    hist(Y1,prob=TRUE,main=varname,xlab="Resid")   
    points(uu,h,type="l")
    
    plot(Y1[-length(Y1)],Y1[-1],pch=3,col="blue",xlim=c(-3.3,3.3),ylim=c(-3.3,3.3),
         xlab=paste(varname,"residual"),ylab=paste(varname,"residual"))
    contour(uu,uu,g,levels=c(0.0025,0.005,0.01,0.02,0.04,0.08,0.12,0.24,0.48),add=TRUE,col="red")
    
    dev.off()
  }
  if(sameday){
    V1  = dimens[1]
    V2  = dimens[2]
    varname1 = names(wacsdata$data)[4+V1]
    varname2 = names(wacsdata$data)[4+V2]
    h1   = wacsdensity_1d(V1,uu,wt,Season.par,Prob,vdry)
    h2   = wacsdensity_1d(V2,uu,wt,Season.par,Prob,vdry)
    g    = wacsdensity_2d(V1,V2,uu,wt,Season.par,Prob,vdry)
    
    Y1 = Y[B,V1] 
    Y2 = Y[B,V2]
    
    if (dry){
      pdf(paste(DIR,"CSN_Season_",s,"_dry_",varname1,"_",varname2,".pdf",sep=""))
    }else{
      pdf(paste(DIR,"CSN_Season_",s,"_wet_",varname1,"_",varname2,".pdf",sep=""))
    }
    par(mfrow=c(2,2))
    
    hist(Y1,prob=TRUE,main=varname1,xlab="Resid")   
    points(uu,h1,type="l")
    
    plot(Y1,Y2,pch=3,col="blue",xlim=c(-3.3,3.3),ylim=c(-3.3,3.3),xlab=paste(varname1,"residual"),ylab=paste(varname2,"residual"))
    contour(uu,uu,g,levels=c(0.0025,0.005,0.01,0.02,0.04,0.08,0.12,0.24,0.48),add=TRUE,col="red")
    
    plot(Y2,Y1,pch=3,col="blue",xlim=c(-3.3,3.3),ylim=c(-3.3,3.3),xlab=paste(varname2,"residual"),ylab=paste(varname1,"residual"))
    contour(uu,uu,t(g),levels=c(0.0025,0.005,0.01,0.02,0.04,0.08,0.12,0.24,0.48),add=TRUE,col="red")
    
    hist(Y2,prob=TRUE,main=varname2,xlab="Resid")   
    points(uu,h2,type="l")
    dev.off()

  }
}
   



