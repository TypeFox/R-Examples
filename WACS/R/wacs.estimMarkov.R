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

  wacs.estimMarkov=function(prob,clustering="hard"){
  ###################################################
  #
  # wacs.estimMarkov: Estimates the 1st order Markov transitions between 
  #                   clusters for one given season
  #
  # ARGUMENTS
  #    prob       : probabilities of belonging to each weather cluster
  #    clustering : indicates wheter 'hard' or 'soft' classficiation must be considered 
  #
  # VALUE
  #    A matrix of transition   
  #
  ###################################################

  if ( (clustering != "hard") & (clustering != "soft") ) stop("[estim.markov] clustering must be 'hard' or 'soft' ")
  
  nDays         = dim(prob)[1]
  nWeatherState = dim(prob)[2]
  
  if (clustering == "hard"){
    WeatherState = apply(prob,1,which.max)
    prob = matrix(0,nDays,nWeatherState)
    for (j in 1:nDays){ 
      prob[j,WeatherState[j]] = 1
    }
  }
  

  M = matrix(0,nWeatherState,nWeatherState)
  if (clustering=="hard") {
    for (j in 1:(nDays-1)){
      M[WeatherState[j],WeatherState[j+1]]=M[WeatherState[j],WeatherState[j+1]] + 1
    }  
    M=M/apply(M,1,sum)
  }else{}
  
  if (clustering=="soft") {
    for (j in 1:(nDays-1)){
      M = M + prob[j,] %o% prob[j+1,]
    }  
    M=M/apply(M,1,sum)
  }else{}
  return(M)
}


