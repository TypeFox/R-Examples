
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

dde.syst <- function(t=NULL,
                     y=NULL,
                     ic=NULL,
                     parms=NULL,
                     derparms=NULL,
                     delags=NULL,
                     codedde=NULL,
                     dosing=NULL,
                     has.dosing=NULL,
                     dose.states=NULL,
                     xdata=NULL,
                     covdata=NULL,
                     scale=NULL,
                     check=FALSE){

  # Input validation
  if (check){
    if (is.null(codedde))
      stop('codedde argument is NULL.')
    
    if (!is.character(codedde))
      stop('codedde argument must be an object of class \'character\'')
    
    if (length(codedde)>1)
      stop('codedde argument must have a length of 1.')
  }
  
  # Evaluates system at past times
  lags <- delags
  
  names(lags) <- paste('alag',names(delags),sep='.')
  t0 <- xdata[1]
  ic <- ic
  
  ylag <- lapply(lags,
    function(x,...){
      if (t-x>=t0) {
        lagvalue(t-x)
      } else {
        ic
      }
    },
    t,t0,ic)
  rm(lags,t0,ic)
  
  # Evaluation of codedde
  cparms <- c(y,ylag,parms,derparms,delags,
              covdata)
  
  dadt <- with(cparms,{
    
    eval(parse(text=codedde))
    
    return(dadt)
    
  })
 
  if (check){
    if (is.null(dim(dadt)[1]))
      stop('dadt in $DDE is not a matrix of dimension (1 x s).')
  }
 
  # Get the variable size info and does some comparisons
  nstate <- dim(dadt)[1]
  
  # Initialize input
  input <- rep(0,nstate)
  
  # Build input
  if (has.dosing){
    for (i in dose.states){
      stdosing <- dosing[dosing[,2]==i,]
      input[i] <- approx(x=stdosing[,1],
                         y=stdosing[,4],
                         xout=t,
                         method='constant',
                         yleft=0,
                         yright=stdosing[dim(stdosing)[1],4],
                         ties='ordered')$y
    }
  }
  
  if (check){
    if (dim(dadt)[1]!=length(scale))
      stop('dadt must have the same dimensions as init and scale.')
  }
  
  # Add the input to the ode system
  dadt <- dadt + input/scale
  
  return(list(c(dadt)))
  
}