
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

ode.model <- function(parms=NULL,
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
    xdata <- xdata
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
  
  # Update initial conditions with bolus dosing if necessary
  sol <- do.call(ode,
                 c(list(y=ic,
                        times=tspan[,1],
                        func=ode.syst,
                        parms=parms,
                        derparms=derparms,
                        codeode=code$de,
                        dosing=dosing,
                        has.dosing=has.dosing,
                        dose.states=dose.states,
                        covdata=covdata,
                        scale=scale,
                        check=check),
                   method=method,
                   options))
  
  ic  <- init.update(a=sol[sol[,1]==tspan[1,1],],
                     t=tspan[1,1],
                     dosing=dosing,
                     scale=scale)
  
  # Integration
  f <- NULL
  for (i in 1:nintervals) {
    # Evaluation times
    eval.times <- xdata[xdata>=tspan[1,i] & xdata<=tspan[2,i]]
    if (is.element(tspan[1,i],eval.times)){
      is.mintspan.in.xdata <- TRUE
    } else {
      is.mintspan.in.xdata <- FALSE
      eval.times <- c(tspan[1,i],eval.times)
    }
    if (is.element(tspan[2,i],eval.times)){
      is.maxtspan.in.xdata <- TRUE
    } else {
      is.maxtspan.in.xdata <- FALSE
      eval.times <- c(eval.times,tspan[2,i])
    }
    
    # Evaluate the solution within the intervals and assumes no observation at bolus times
    sol <- do.call(ode,
                   c(list(y=ic,
                          times=eval.times,
                          func=ode.syst,
                          parms=parms,
                          derparms=derparms,
                          codeode=code$de,
                          dosing=dosing,
                          has.dosing=has.dosing,
                          dose.states=dose.states,
                          covdata=covdata,
                          scale=scale,
                          check=FALSE),
                     method=method,
                     options))
    
    # initialize states for next loop iteration
    if (i!=nintervals){
      ic <- init.update(a=sol[sol[,1]==tspan[2,i],],
                        t=tspan[1,i+1],
                        dosing=dosing,
                        scale=scale)
    }
    
    # Filter sol based upon is.mintspan.in.xdata and is.maxtspan.in.xdata
    if (!is.mintspan.in.xdata){
      sol <- sol[-1,]
    } else {
      if (!is.null(f)) f <- f[-size(f,1),]
    }
    if (!is.maxtspan.in.xdata){
      ftmp <- sol[-size(sol,1),]
    } else {
      ftmp <- sol
    }
    
    # Concatenate Ftmp to the previous predictions
    if (is.null(f)) {
      f <- ftmp
    } else {
      f <- rbind(f,ftmp)
    }
  }
  
  # Define ouput from the system
  f <- de.output(f=transpose(f),
                 parms=parms,
                 derparms=derparms,
                 codeoutput=code$output,
                 dosing=dosing,
                 xdata=xdata)
  
  # Re-attach evaluation times xdata for simulation run only
  if (issim > 0.5){
    f <- rbind(xdata,f)
  }
  
  return(f)
  
}
