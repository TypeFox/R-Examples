
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

explicit.model <- function(parms=NULL,
                           derparms=NULL,
                           code=NULL,
                           bolus=NULL,
                           infusion=NULL,
                           xdata=NULL,
                           covdata=NULL,
                           issim=0,
                           check=FALSE){
  
  # Input validation
  if (check){
    if (length(parms)==0)
      stop('parms argument is NULL.')
    
    if (size(bolus,2)!=4)
      stop('bolus argument does not have a di x 4 dimesion.')
    
    if (size(infusion,2)<4)
      stop('infusion argument contain less than 4 columns.')
    
    if (length(xdata)<2)
      stop('xdata argument must contain at least 2 elements.')
  }
  
  # Sort dosing by time
  bolus <- bolus[order(bolus[,1]),]
  infusion <- infusion[order(infusion[,1]),]
  
  # Determine integration intervals
  tspan <- create.intervals(xdata=xdata,
                            bolus=bolus,
                            infusion=infusion)
  nintervals <- size(tspan,2)
  
  # Determine the time points for model evaluation
  if (issim < 0.5){
    xdata <- xdata
  } else {
    xdata <- NULL
    nint <- ceiling(1001/nintervals)
    # Check that nint is odd; if not, adds 1
    if (!nint%%2)
      nint <- nint + 1
    
    # Create vector of time
    for (i in 1:nintervals){
      xtmp <- seq(tspan[1,i],tspan[2,i],length.out=nint)
      if (i==1){
        xdata <- c(xdata,xtmp)
      } else {
        xdata <- c(xdata,xtmp[2:length(xtmp)])
      }
    }
  }
  
  # Create time variable from xdata
  times <- xdata
  
  # Evaluation of codeic
  if (is.null(covdata)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,covdata)
  }
  
  f <- with(cparms,{
    
    eval(parse(text=code$output))
    
    return(y)
    
  })
  
  # Re-attach evaluation times for simulation run only
  if (issim > 0.5){
    f <- rbind(times,f)
  }
  
  return(f)
  
}

