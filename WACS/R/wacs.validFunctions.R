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
# but WITHOUT ANY WARRANTY; without even the implied warrSanty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################

#
#  Specific Functions for the validation of WACS simulations, compared to WACS data
#

##################### One Simulation #####################

wacsvalid.Sim = function(wacsdata,wacspar,wacssimul,varname){
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid.Sim] Warning: Nb days simulated should be 
                a multiple of 365")
  }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")
  
  # All usual variable, including 'tmax'
  if (varname !="tmoy"){
    y = extract.annual.trend(sims[,varname],spar=wacspar$Trend$Param[1],trend.norm=wacspar$Trend$Param[2])
    CentralSim   = y[,1]
    DeviationSim = y[,2]
    
    if (wacsdata$Trange && (varname=="tmax")){
      Obs = wacsdata$data[,"tmin"] + wacsdata$data[,"trange"]
      y = extract.annual.trend(Obs,spar=wacspar$Trend$Param[1],trend.norm=wacspar$Trend$Param[2])
      CentralObs   = y[,1]
      DeviationObs = y[,2]
    }else{
      Obs = wacsdata$data[,varname]
      y = extract.annual.trend(Obs,spar=wacspar$Trend$Param[1],trend.norm=wacspar$Trend$Param[2])
      CentralObs   = y[,1] 
      DeviationObs = y[,2] 
    }
    
    zObs = matrix(0,NyObs,365)
    zSim = matrix(0,NyObs,365)
    for (i in 1:NyObs){
      zObs[i,] = wacsdata$data[((i-1)*365+1):(i*365),varname]
    }
    for (i in 1:NySim){
      zSim[i,] = sims[((i-1)*365+1):(i*365),varname]
    }
  }else{   # If the variable is 'tmoy'
    
    zObs = matrix(0,NyObs,365)
    zSim = matrix(0,NyObs,365)
    if (wacsdata$Trange){
      Obs = wacsdata$data[,"tmin"] + wacsdata$data[,"trange"]/2
      Sim = sims[,"tmin"] + sims[,"trange"]/2
    }else{
      Obs = (wacsdata$data[,"tmin"] + wacsdata$data[,"tmax"])/2
      Sim = (sims[,"tmin"] + sims[,"tmax"])/2
    }
    for (i in 1:NyObs){      
      zObs[i,] = Obs[((i-1)*365+1):(i*365)]
    }
    for (i in 1:NySim){
      zSim[i,] = Sim[((i-1)*365+1):(i*365)]
    }
    
    y = extract.annual.trend(Obs,spar=wacspar$Trend$Param[1],trend.norm=wacspar$Trend$Param[2])
    CentralObs   = y[,1]
    DeviationObs = y[,2]
    y = extract.annual.trend(Sim,spar=wacspar$Trend$Param[1],trend.norm=wacspar$Trend$Param[2])
    CentralSim   = y[,1]
    DeviationSim = y[,2]
  }
  
  res=list(varname=varname,CentralObs=CentralObs,DeviationObs=DeviationObs,
           zObs=zObs,CentralSim=CentralSim,DeviationSim=DeviationSim,zSim=zSim,seasons=wacsdata$seasons)
  class(res) = "WACSvalidSim"
  return(res)
  
}

##################### Rain #####################

wacsvalid.Rain = function(wacsdata,wacspar){
  nbSeasons = length(wacspar$seasons$day);
  res = list();
  for (s in 1:nbSeasons){
    scale = wacspar$Rain$RainPar[s,1];
    shape = wacspar$Rain$RainPar[s,2];
    y     = sort(wacsdata$data$rain[(wacsdata$data$season == s) & 
                                (wacsdata$data$rain > 0)]);
    res[[s]] = list(
      theoretical = qgamma(c(1:(length(y)-1))/length(y),
                         scale=scale, shape=shape),  
      observed = y[1:(length(y)-1)], par=wacspar$Rain$RainPar[s,]);
  }
  class(res) = "WACSvalidRain";
  return(res)
}

##################### MeanSd #####################

wacsvalid.MeanSd = function(wacsdata,wacssimul,varname){
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | 
                               (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid] Warning: Nb days simulated should be 
                a multiple of 365")
  }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")
  
  meanObs  =  matrix(0,NyObs,12)
  sdObs    =  matrix(0,NyObs,12)
  meanSim  =  matrix(0,NySim,12)
  sdSim    =  matrix(0,NySim,12)
  
  for (i in 1:NyObs) {
    y = unique(wacsdata$data$year)[i]
    for ( j in 1:12 ) {
      meanObs[i,j] = mean(wacsdata$data[(wacsdata$data$year == y 
                                         & wacsdata$data$month == j), varname])
      sdObs[i,j] = sd(wacsdata$data[(wacsdata$data$year == y 
                                     & wacsdata$data$month == j), varname])
    }
  }

  for ( i in 1:NySim) {
    y = unique(sims$year)[i]
    for ( j in 1:12 ) {
      meanSim[i,j] = mean(wacssimul$sim[(wacssimul$sim$year == y 
                                         & wacssimul$sim$month == j), varname]);
      sdSim[i,j] = sd(wacssimul$sim[(wacssimul$sim$year == y 
                                     & wacssimul$sim$month == j), varname]);
    }
  }
  res = list(meanObs=meanObs, sdObs=sdObs, meanSim=meanSim, sdSim=sdSim,
             varname=varname)
  class(res) = "WACSvalidMeanSd";
  return(res)
}

##################### BiVar #####################

wacsvalid.BiVar = function(wacsdata,wacssimul,varname,varname2){
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | 
                               (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid.BiVar] Nb days simulated should be 
          a multiple of 365")
    }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")
  
  month    = 1:12
  
  corObs   = matrix(0,NyObs,12)
  corSim   = matrix(0,NySim,12)

  for ( i in 1:NyObs) {
    y = sort(unique(wacsdata$data$year))[i]
    for ( j in 1:12) {
        tmp  = wacsdata$data[which((wacsdata$data$year == y) &
                                    (wacsdata$data$month == j)),varname]
        tmp2 = wacsdata$data[which((wacsdata$data$year == y) &
                                     (wacsdata$data$month == j)),varname2]
        corObs[i,j] = cor(tmp, tmp2)
      }
    }
  
    for ( i in 1:NySim) {
      y = sort(unique(sims$year))[i]
      for ( j in 1:12) {
        tmp  = sims[which((sims$year == y) & (sims$month == j)), varname]
        tmp2 = sims[which((sims$year == y) & (sims$month == j)), varname2]
        corSim[i,j] = cor(tmp, tmp2)
      }
    }
    res = list(corObs= corObs, corSim=corSim, varname=varname,varname2=varname2);
    class(res) = "WACSvalidBiVar";
    return(res)
}

##################### Temporal Correlation #####################

wacsvalid.CorTemp = function(wacsdata,wacssimul,varname){
  options(warn=-1)
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | 
                               (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid.CorTemp] for 'CorTemp' nb days simulated should be 
                          a multiple of 365")
  }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")

  month    = 1:12
  
  corObs =  matrix(0,NyObs,12)
  corSim =  matrix(0,NySim,12)
  
  for ( i in 1:NyObs) {
    y = unique(wacsdata$data$year)[i]

    for ( j in 1:12) {
      tmp = wacsdata$data[which((wacsdata$data$year == y) & (wacsdata$data$month == j)), varname]
      corObs[i,j] = cor(tmp[-1], tmp[-length(tmp)])
    }
  }

  for ( i in 1:NySim) {
    y = unique(sims$year)[i]
    for ( j in 1:12) {
      tmp = sims[which((sims$year == y) & (sims$month == j)), varname]
      corSim[i,j] = cor(tmp[-1], tmp[-length(tmp)])
    }
  }
  options(warn=0)
  res = list(corObs= corObs, corSim=corSim, varname=varname);
  class(res) = "WACSvalidCorTemp";
  return(res)
}

##################### SumBase #####################

wacsvalid.SumBase = function(wacsdata,wacssimul,varname,base=0,months=1:12){
  
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | 
                               (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid] Nb days simulated should be 
                a multiple of 365")
  }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")
  
  SumObs =  rep(0,NyObs)
  SumSim =  rep(0,NySim)
  
  if (varname=="tmoy"){
    tmoyObs = (wacsdata$data$tmin + wacsdata$data$tmax)/2
    tmoySim = (sims$tmin  + sims$tmax)/2
    for (i in 1:NyObs) {
      y      = unique(wacsdata$data$year)[i]
      selObs = (wacsdata$data$year == y) & (wacsdata$data$month %in% months) & (tmoyObs > base)
      SumObs[i] = sum(tmoyObs[selObs]) 
    }  
    for (i in 1:NySim) {
      y      = unique(sims$year)[i]
      selSim = (sims$year == y) & (sims$month %in% months) & (tmoySim > base)
      SumSim[i] = sum(tmoySim[selSim])
    }
  }else{
    for (i in 1:NyObs) {
      y      = unique(wacsdata$data$year)[i]
      selObs = (wacsdata$data$year == y) & (wacsdata$data$month %in% months) & (wacsdata$data[,varname] > base)
      SumObs[i] = sum(wacsdata$data[selObs,varname]) 
    }  
    
    for (i in 1:NySim) {
      y      = unique(sims$year)[i]
      selSim = (sims$year == y) & (sims$month %in% months) & (sims[,varname] > base)
      SumSim[i] = sum(sims[selSim,varname])
    }
  }
  
  res = list(SumObs=SumObs, SumSim=SumSim, varname=varname,base=base)
  class(res) = "WACSvalidSumBase"
  return(res)
}

##################### Persistence

wacsvalid.Persistence = function(wacsdata,wacssimul,varname,base=0,above=TRUE,months=1:12){
  
  sims = wacssimul$sim[which((wacssimul$sim$month != 2) | 
                               (wacssimul$sim$day != 29)),]
  if (nrow(sims) %% 365 != 0) {
    stop ("[wacsvalid.Persistence] Nb days simulated should be 
                a multiple of 365")
  }
  if (varname == "tmoy"){
    stop ("[wacsvalid.Persistence] 'tmoy' cannot be chosen as a variable")
  }
  
  NyObs       = length(unique(wacsdata$data$year))
  NySim       = length(unique(sims$year))
  if (NyObs != NySim) stop("[wacsvalid]  Data and Simulations have different length")
  
  years = unique(wacsdata$data$year)
  z  = wacsdata$data[(wacsdata$data$year == years[1]) & (wacsdata$data$month %in% months),varname]
  varObs = z
  for (i in 2:NyObs) {
    z = wacsdata$data[(wacsdata$data$year == years[i]) & (wacsdata$data$month %in% months),varname]
    varObs = rbind(varObs,z)
  }  

  years = unique(sims$year)
  z  = sims[(sims$year == years[1] ) & (sims$month %in%months), varname]
  varSim = z
  for (i in 2:NySim){
    z  = sims[(sims$year == years[i] ) & (sims$month %in%months), varname]
    varSim = rbind(varSim,z)
  }
  
  FreqObs = persistence(varObs,base,above,months)
  FreqSim = persistence(varSim,base,above,months)
  Persmax = max(which(FreqObs>0),which(FreqSim>0))
  FreqObs = FreqObs[1:Persmax]
  FreqSim = FreqSim[1:Persmax]
 
  res = list(FreqObs=FreqObs, FreqSim=FreqSim, varname=varname, base=base, above=above)
  class(res) = "WACSvalidPersistence"
  return(res)
}

persistence = function(Var,base,above,months) {
  
  #####################################################
  #
  # WACSgen project v2013. Author D. Allard
  #
  #   Function persistence : internal function to compare the persistence of a variable above or below a base
  #
  #   ARGUMENTS :
  #    Var   :         variable to be analyzed; it is an array; each line is a separate year
  #    base  :         threshold
  #    above :         persistence above threshold if TRUE; below threshold if FALSE
  #    months:         Months to be considered
  #
  #
  
  Ny        = dim(Var)[1]
  MaxLength = dim(Var)[2]
  Freq      = rep(0,MaxLength)
  
  for (y in 1:Ny){
    persvar = rep(1,MaxLength)
    for (i in 2:MaxLength){
      if (above){
        if ((Var[y,i] > base) && (Var[y,i-1] > base) ){
          persvar[i] = persvar[i-1] + 1 
          if (i==MaxLength){
            Freq[persvar[i]] = Freq[persvar[i]] + 1
          }
        }
        if ((Var[y,i] <= base) && (Var[y,i-1] > base) ){
        Freq[persvar[i-1]] = Freq[persvar[i-1]] + 1
        }
      }else{
        if ((Var[y,i] <= base) && (Var[y,i-1] <= base) ){
          persvar[i] = persvar[i-1] + 1 
          if (i==MaxLength){
            Freq[persvar[i]] = Freq[persvar[i]] + 1
          }
        }    
        if ((Var[y,i] > base) && (Var[y,i-1] <= base) ){
          Freq[persvar[i-1]] = Freq[persvar[i-1]] + 1
        }
      }
    }
  }
  return(Freq)
}


