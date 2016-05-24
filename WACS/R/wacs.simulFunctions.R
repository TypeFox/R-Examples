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


# Contains functions for performing simulations based on estimated parameters of the WACS model

rMarkov = function(state,TransM){
  ###################################################################
  #
  # rMarkov : simulates a new state according to a transition probability
  # matrix
  #
  # ARGUMENTS
  #      state      : current state
  #      TransM     : transition probability matrix
  #
  # VALUE
  #     newstate : the new state
  ###################################################################
  
  if (state > dim(TransM)[1]) stop ("[rMarkov] state > dim(TransM)")
  prob = TransM[state,]
  newstate = sample(1:length(prob),1,prob=prob)
  return(newstate)
}

map.wt=function(X,WACSparam.s){
  ###############################################
  #
  # map.wt      : finds the most likely weather state for a given vector of 
  #               residuals
  #
  # ARGUMENTS
  #         X     : the vector of residuals
  #   WACSparam.s : parameters for a given season, as created by WACSestim 
  #
  # VALUE
  #            wt : most likely weather type
  #
  # Calls 
  #     dmcsnstar 
  #
  # New version 8 april 2014 
  ###############################################"
  
  # Number of weather state in season s
  Nwt = dim(WACSparam.s$TransM)[1]
  
  # computes the density associated to each weather state in season S;
  # If rain status is not similar, density is set to -999
  
  d      = rep(-999,Nwt)
  WET    = (X[1] > -999)
  if (!WET) {
    wt.ini = 1
    wt.end = WACSparam.s$NumbWT[1]
  }
  if (WET){
    wt.ini = WACSparam.s$NumbWT[1] + 1
    wt.end = WACSparam.s$NumbWT[1] + WACSparam.s$NumbWT[2]
  }
  p.st   = abs(eigen(t(WACSparam.s[[3]]))$vectors[,1]) # finds the stationary distribution of the Markov Chain
  p.st   = p.st/sum(p.st)
  for (w in wt.ini : wt.end){
    Mu    =  WACSparam.s[[w+3]]$loc     			
    Sigma =  WACSparam.s[[w+3]]$cov		
    Skew  =  WACSparam.s[[w+3]]$skew
    for (i in 1:length(Skew)) Skew[i]=sign(Skew[i])*min(0.99,abs(Skew[i]))
    if(!WET) d[w] = dmcsnstar(X[-1],Mu,Sigma,Skew)*p.st[w]
    if(WET) d[w] = dmcsnstar(X,Mu,Sigma,Skew)*p.st[w]
  }
  wt = which(d==max(d))[1]
  #  if (length(wt) > 1) {
  #      wt = wt[1]
  #  }
  return(wt)
}

boundsvariables = function(X,bounds,iitmin,iitrange){
  ###################################################################
  #
  # Postprocess the variables to force them within accepted bounds
  #
  # ARGUMENT 
  #           X: vector of variables
  #     vbounds: lower and upper bounds for each variables
  #     Trange: = TRUE   : X[1] is Tmin; X[2] is Tmax - Tmin
  #              = FALSE  : X[1] is Tmin; X[2] is Tmax  ==> condition is to check is Tmax > Tmin
  #
  # VALUE
  #          Res: new vector of variables
  ###################################################################
  
  Nv = length(X)
  Res = X
  for (v in 1:Nv){
    if (X[v] < bounds[1,v]) Res[v] = bounds[1,v]
    if (X[v] > bounds[2,v]) Res[v] = bounds[2,v]
  }
  if (iitrange != -1) {
    #use extra bound for Trange
    if ((X[iitmin] + X[iitrange]) < bounds$tmax[1]) {
      Res[iitrange] = bounds$tmax[1] - X[iitmin]
    }
    if ((X[iitmin] + X[iitrange]) > bounds$tmax[2]) {
      Res[iitrange] = bounds$tmax[2] - X[iitmin]
    }
  }
  if (iitrange == -1){
    if (X[iitmin] > X[iitmin+1]){
      Res[iitmin+1] = Res[iitmin]
    }
  }
  return (Res);  
}

testvariables = function(X,bounds,iitmin,iitrange){
  ###################################################################
  #
  # Checks wether the variables are within accepted bounds
  #
  # ARGUMENT 
  #           X: vector of variables
  #     vbounds: lower and upper bounds for each variables
  #     Trange: = TRUE   : X[1] is Tmin; X[2] is Tmax - Tmin
  #              = FALSE  : X[1] is Tmin; X[2] is Tmax  ==> condition is to check is Tmax > Tmin
  #
  # VALUE
  #          OK: a Boolean
  ###################################################################
  ACCEPT = TRUE
  Nv = length(X)
  for (v in 1:Nv){
    if ((X[v] < bounds[1,v]) || X[v] > bounds[2,v]) {    
      ACCEPT = FALSE
    }
  }
  if (iitrange != -1) {
    #use extra bound for Trange
    if ((X[iitmin] + X[iitrange]) < bounds$tmax[1]) {
      ACCEPT = FALSE;
    }
    if ((X[iitmin] + X[iitrange]) > bounds$tmax[2]) {
      ACCEPT = FALSE;
    }
  }
  if (iitrange == -1){
    if (X[iitmin] > X[iitmin+1]){
      ACCEPT = FALSE
    }
  }
  return (ACCEPT);  
}


transform.rain = function(y,rain.model,rain.par){
  ###################################################################
  #
  #
  # Transforms Gaussian scores into rain values
  #
  # ARGUMENT 
  #           y: Gaussian score
  #  rain.model: "Gamma" for a Gamma distribution; "A-B" for the Allard-Bourotte model 
  #  rain.par  : parameters of the model
  # VALUE
  #           z: a transformed value
  ###################################################################
  if(rain.model=="Gamma"){
    z = qgamma(pnorm(y),scale = rain.par[1],shape=rain.par[2])
  } else if (rain.model=="AB"){
    stop ("[WACSsimul] transform.rain: model AB not implemented yet")
  }
  if(rain.model=="None"){
    z = y
  }
  return(z)  
}


