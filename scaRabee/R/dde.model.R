
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

dde.model <- function(parms=NULL,
                      derparms=NULL,
                      code=NULL,
                      bolus=NULL,
                      infusion=NULL,
                      xdata=NULL,
                      covdata=NULL,
                      issim=0,
                      check=FALSE,
                      options=list(method='lsoda')){
  
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
  
  # Process dosing information
  dosing <- make.dosing(allparms=c(derparms,as.list(parms)),
                        bolus=bolus,
                        infusion=infusion,
                        check=check)
  
  has.dosing <- any(dosing[,4]>0)
  dose.states <- unique(dosing[dosing[,4]>0,2])
  
  # Determine integration intervals
  tspan <- create.intervals(xdata=xdata,
                            bolus=bolus,
                            infusion=infusion)
  nintervals <- size(tspan,2)
  
  # Determine the time points for model evaluation
  if (issim < 0.5){
    xdata.ori <- xdata
    if (tspan[1,1]<xdata[1]){
      xdata <- c(tspan[1,1],xdata)
    }
    
  } else {
    obst <- xdata
    xdata <- NULL
    nint <- ceiling(1001/nintervals)
    # Checks that nint is odd; if not, adds 1
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
    xdata <- sort(unique(c(xdata, obst)))
    xdata.ori <- xdata
  }
  
  # Extract method, and other options from options
  method <- options$method
  options$method <- NULL
  
  # Define initial conditions
  ic <- init.cond(parms=parms,
                  derparms=derparms,
                  codeic=code$ic,
                  dosing=bolus,
                  check=check)
  
  # Determine the scaling factors for inputs
  scale <- input.scaling(parms=parms,
                         derparms=derparms,
                         codescale=code$scale,
                         ic=ic,
                         check=check)
  
  # Determine system events
  events <- get.events(bolus=bolus,
                       scale=scale)
  
  # Get the delay parameters
  lags <- dde.lags(parms=parms,
                   derparms=derparms,
                   codelags=code$lag,
                   check=check)
  
  # Check dde system
  if (check){
    if (dim(events)[1]>0){
      sol <- do.call(dede,
                     c(list(y=ic,
                            times=c(0,0.1),
                            func=dde.syst,
                            ic=ic,
                            parms=parms,
                            derparms=derparms,
                            delags=lags,
                            codedde=code$de,
                            dosing=dosing,
                            has.dosing=has.dosing,
                            dose.states=dose.states,
                            xdata=xdata,
                            covdata=covdata,
                            scale=scale,
                            check=check,
                            events=list(data=events)),
                       method=method,
                       options))
    } else {
      sol <- do.call(dede,
                     c(list(y=ic,
                            times=c(0,0.1),
                            func=dde.syst,
                            parms=parms,
                            ic=ic,
                            derparms=derparms,
                            delags=lags,
                            codedde=code$de,
                            dosing=dosing,
                            has.dosing=has.dosing,
                            dose.states=dose.states,
                            xdata=xdata,
                            covdata=covdata,
                            scale=scale,
                            check=check),
                       method=method,
                       options))
    }
    check <- FALSE
  }
  
  # Solve DDE system (update initial conditions with bolus dosing if necessary
  if (dim(events)[1]>0){
    sol <- do.call(dede,
                   c(list(y=ic,
                          times=xdata,
                          func=dde.syst,
                          ic=ic,
                          parms=parms,
                          derparms=derparms,
                          delags=lags,
                          codedde=code$de,
                          dosing=dosing,
                          has.dosing=has.dosing,
                          dose.states=dose.states,
                          xdata=xdata,
                          covdata=covdata,
                          scale=scale,
                          check=check,
                          events=list(data=events)),
                     method=method,
                     options))
  } else {
    sol <- do.call(dede,
                   c(list(y=ic,
                          times=xdata,
                          func=dde.syst,
                          ic=ic,
                          parms=parms,
                          derparms=derparms,
                          delags=lags,
                          codedde=code$de,
                          dosing=dosing,
                          has.dosing=has.dosing,
                          dose.states=dose.states,
                          xdata=xdata,
                          covdata=covdata,
                          scale=scale,
                          check=check),
                     method=method,
                     options))
  }
  
  # Define ouput from the system
  f <- de.output(f=transpose(as.matrix(sol)),
                 parms=parms,
                 derparms=derparms,
                 codeoutput=code$output,
                 dosing=dosing,
                 xdata=xdata)
  
  # Filter f to extract only time points in xdata.ori
  f <- f[,match(xdata.ori,sol[,1]),drop=FALSE]
  
  # Re-attach evaluation times xdata for simulation run only
  if (issim > 0.5){
    f <- rbind(xdata,f)
  }
  
  return(f)
  
}
