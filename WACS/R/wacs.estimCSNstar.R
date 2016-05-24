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
  
  wacs.estimCSNstar=function(Y,z,jday){
  ###################################################################
  #
  #
  # wacs.estimCSNstar: does the estimation of the parameters of a multivariate
  #                    CSN* Random Variable, according to a weighted moment method. 
  #                    Flecher, Naveau and Allard (2009), Statistics and 
  #                    Probability Letters, 79, 1977-1984.
  #
  # ARGUMENT 
  #            Y: matrix (nxp) of values
  #            z: vector of weigths = probability of belonging to the class
  #         jday: julian jday (allows to know if jdays are consecutive)
  #
  # VALUE
  #    A list containing the location vector, the covariance matrix,
  #    and the shape parameter vector. Shape parameters are between (-1 and 1).
  #    Maximum skewness is obtained with -1 or 1. A value of 0 corresponds to
  #    no skewness
  #
  # CALLS
  #    csnPWM
  #    WMdiff
  #
  # New version 8 april 2014
  ###################################################################
  

  # Number of variables and number of data
  Nv = dim(Y)[2]
  Nd = dim(Y)[1]
  
  # For robust estimation
  for (i in 1:Nd){
    for (v in 1:Nv){
      if (Y[i,v] > 3.5){ Y[i,v] = 3.5}
      if (Y[i,v] < -3.5){ Y[i,v] = -3.5}
    }
  }
  
  # Computation of first, second and third moments
  E = rep(0,Nv)
  for (v in 1:Nv){ 
    E[v] = sum(Y[,v]*z)
  }
  E = E/sum(z)    

  V = matrix(0,ncol=Nv,nrow=Nv)
  for (v in 1 : Nv){
    for (w in 1 : Nv){
      V[v,w] = sum(Y[,v]*Y[,w]*z)/sum(z)
      V[v,w] = V[v,w] - E[v]*E[w]
    }
  }
  
  M0 = 0
  for (i in 1:dim(Y)[1]){ 
    M0 = M0 + pmnorm(Y[i,],varcov=diag(Nv))[1]*z[i]
  }
  M0 = M0/sum(z)
  
  # Estimating the correlation coefficients for two consecutive jdays
  EE = rep(0,Nv); UU = EE; WW  = EE
  Consec.jdays = (jday[-1]-jday[-Nd]==1)
  for (v in 1:Nv){
    EE[v] = sum(Y[-Nd,v]*z[-1]*z[-Nd]*Consec.jdays)/sum(z[-1]*z[-Nd]*Consec.jdays)   
    UU[v] = sum(Y[-Nd,v]^2*z[-Nd]*z[-1]*Consec.jdays)/sum(z[-Nd]*z[-1]*Consec.jdays) 
    WW[v] = sum(Y[-Nd,v]*Y[-1,v]*z[-Nd]*z[-1]*Consec.jdays)/sum(z[-Nd]*z[-1]*Consec.jdays)
    UU[v] = UU[v] - EE[v]^2
    WW[v] = WW[v] - EE[v]^2
  }   
  
  if(is.na(EE[1])) stop ("[WACSestim] not enough jdays in one of the WT for hard classification")
  Rho.hat = WW/UU
  for (v in 1:Nv) {
    if (Rho.hat[v] >  1) Rho.hat[v] = 0.95
    if (Rho.hat[v] < -1) Rho.hat[v] = -0.95
  }
  
  
  # Initial "guesstimate" of Shat
  Shat.ini=c()  
  for (v in 1:Nv){ 
    Shat.ini=c(Shat.ini,csnPWM(Y[,v],z)[[3]])
  }
  
  # Estimating the skewness parameters by profile MOM
  Shat = optim(par=Shat.ini,fn=WMdiff,moments=c(E,V,Rho.hat,M0),method="Nelder-Mead",
               control=list(abstol=0.001))$par
  for (v in 1:Nv){
    if (Shat[v] >= 0.99 )  Shat[v] =  0.9
    if (Shat[v] <= -0.99 ) Shat[v] = -0.9
  } 
  
  
  Sigma.hat = sqrtm(V)%*%solve(diag(Nv) - 2/pi*diag(Shat)^2)%*%sqrtm(V)

  
  # Creating SS and imposing positive definitness
  SS        = matrix(0,2*Nv,2*Nv)
  SS[1:Nv,1:Nv]                     = Sigma.hat
  SS[(Nv+1):(2*Nv),(Nv+1):(2*Nv)]   = Sigma.hat
  SS[1:Nv,(Nv+1):(2*Nv)]            = sqrtm(Sigma.hat)%*%diag(Rho.hat)%*%sqrtm(Sigma.hat)
  SS[(Nv+1):(2*Nv),1:Nv]            = t(SS[1:Nv,(Nv+1):(2*Nv)])
  # SS        = def.pos(SS)
  # Sigma.hat = SS[1:Nv,1:Nv]
  # for (i in 1:Nv){ 
  #    Rho.hat[i] = SS[i,(Nv+i)]/SS[i,i]
  # }
  # Estimating the expectations
  Mu2    = rep(E,2) - sqrtm(SS)%*%diag(rep(Shat,2))%*%rep(sqrt(2/pi),2*Nv)
  Mu.hat = 0.5*(Mu2[1:Nv]+Mu2[(Nv+1):(2*Nv)])
  # Mu.hat = E - diag(Shat)%*%sqrtm(Sigma.hat)%*%rep(sqrt(2/pi),Nv)
  
  par.csn=list(loc=as.vector(Mu.hat),cov=Sigma.hat,skew=Shat,rho=Rho.hat)
  return(par.csn)
}
  
#####################################  End of function ##############


WMdiff = function(Shat,moments){
  ###################################################################
  #
  # WMdiff auxiliary function to be minimized
  #
  # ARGUMENT 
  #     Shat   : Skewness values
  #     moments: experimental moments
  #
  # VALUE
  #    diff    : a criterion to be minimized
  #
  ###################################################################
  
  Nv = length(Shat)
  E  = moments[1:Nv]
  V  = matrix(moments[(Nv+1):(Nv*(Nv+1))],ncol=Nv)
  Rho.hat = moments[((Nv*(Nv+1))+1):(Nv*(Nv+2))]
  M0 = moments[(Nv*(Nv+2))+1]
  Sigma.hat = sqrtm(V)%*%solve(diag(Nv)-2/pi*diag(Shat)^2)%*%sqrtm(V)
  

  
  # Imposing positive definitness to SS 
  # 07/04/2014
  # Not necessary anymore in the new formulation of time correlations: 
  # Sigma.hat is d.p. and will always lead to a d.p. matrix SS 
  #
  SS = matrix(0,2*Nv,2*Nv)
  SS[1:Nv,1:Nv]                     = Sigma.hat
  SS[(Nv+1):(2*Nv),(Nv+1):(2*Nv)]   = Sigma.hat
  SS[1:Nv,(Nv+1):(2*Nv)]            = sqrtm(Sigma.hat)%*%diag(Rho.hat)%*%sqrtm(Sigma.hat)
  SS[(Nv+1):(2*Nv),1:Nv]            = t(SS[1:Nv,(Nv+1):(2*Nv)])
  # SS = def.pos(SS)
  # Sigma.hat = SS[1:Nv,1:Nv]
  Mu2    = rep(E,2) - sqrtm(SS)%*%diag(rep(Shat,2))%*%rep(sqrt(2/pi),2*Nv)
  Mu.hat = 0.5*(Mu2[1:Nv]+Mu2[(Nv+1):(2*Nv)])
  
  
  # Estimating the expectation
  # Mu.hat = E - diag(Shat)%*%sqrtm(Sigma.hat)%*%rep(sqrt(2/pi),Nv)

  MM     = diag(2*Nv)
  MM[1:Nv,1:Nv]          = MM[1:Nv,1:Nv] + Sigma.hat
  MM[(Nv+1):(2*Nv),1:Nv] = diag(Shat)%*%sqrtm(Sigma.hat)
  MM[1:Nv,(Nv+1):(2*Nv)] = t(MM[(Nv+1):(2*Nv),1:Nv])
  E0  = 2^Nv * pmnorm(rep(0,2*Nv),mean=c(Mu.hat,rep(0,Nv)),varcov=MM)[1]
  pen = rep(0,Nv)
  for (v in 1:Nv){ 
    pen[v] = 5*(abs(Shat[v])-1)
  }
  diff = abs(E0-M0) + sum(pen)
}

csnPWM <-function(yy,zz)
{
  ###################################################################
  #
  # csnPWM: estimates the shape parameter of a CSN* univariate Random Variable
  #
  # ARGUMENT 
  #     yy: vector of values
  #     zz: weigths
  #
  # VALUE
  #    estimates: a vector of estimates for the three parameters:
  #    location, shape, skewness
  ###################################################################
  
  if (length(yy) != length(zz)){
    stop("[csnPWM] vector of weight and vector of data must have same length")
  }
  
  # computing the moments (order 1 et 2 and weighted moment)
  yy = yy[zz>0.5]
  zz = zz[zz>0.5]
  S   = sum(zz)
  m1  = sum(yy*zz)/S
  m2  = sum(yy^2*zz)/S - m1*m1
  m11 = sum(yy*pnorm(yy)*zz)/S


  output = optim(c(m1,sqrt(m2),0),flik,gr=NULL,yy=yy,zz=zz,method="Nelder-Mead")
  if (output$par[3]> 0.98) {output$par[3] = 0.98}
  if (output$par[3]< -0.98) {output$par[3] = -0.98}
  if (output$par[2]< 0){
    ouptut$par[1] = m1
    output$par[2] = sqrt(m2)
    ouptut$par[3] = 0
  }
            
  estimates = round(output$par,5)
  names(estimates) = c("location","scale","skewness")
  return(estimates)
}

  flik = function(param,yy,zz){
    LL = zz
    for (i in 1:length(yy)) LL[i] = log(zz[i]*dcsn(yy[i], location=param[1],scale=param[2],skew=param[3]))
    LLS = -sum(LL)
    return(LLS)
  }
  



