
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

problem.eval <- function(subproblem=NULL,x=NULL,grid=FALSE,check=FALSE){
  
  # Copy subproblem$init in newparam and replace the new estimates in newparam
  newparam <- subproblem$init
  
  if (grid) {
    index <- which(newparam$isfix==0)
    newparam$value[index] <- x[index]
  } else {
    npar     <- length(newparam$names)
    estindex <- 1
    for (i in 1:npar){
      if (newparam$isfix[i]==0){
        newparam$value[i] <- x[estindex]
        estindex <- estindex + 1
      }
    }
  }
  
  # Retrieve primary parameters
  parms <- c(get.parms.data(x=newparam,which='value',type='P'),
             get.parms.data(x=newparam,which='value',type='L'),
             get.parms.data(x=newparam,which='value',type='IC'),
             get.parms.data(x=newparam,which='value',type='V'))
  names(parms) <- c(get.parms.data(x=newparam,which='names',type='P'),
                    get.parms.data(x=newparam,which='names',type='L'),
                    get.parms.data(x=newparam,which='names',type='IC'),
                    get.parms.data(x=newparam,which='names',type='V'))
  attr(parms,'type') <- c(get.parms.data(x=newparam,which='type',type='P'),
                          get.parms.data(x=newparam,which='type',type='L'),
                          get.parms.data(x=newparam,which='type',type='IC'),
                          get.parms.data(x=newparam,which='type',type='V'))
  
  # Retrieve derived parameters
  if (subproblem$modfun%in%c('ode.model','dde.model')){
    derparms <- derived.parms(parms=parms,
                              covdata=subproblem$cov,
                              codederiv=subproblem$code$deriv,
                              check=check)
  } else {
    derparms <- NULL
  }
  
  # Evaluate the full model function modfun
  if (subproblem$modfun%in%c('ode.model','dde.model')){
    fpred <- do.call(eval(parse(text=subproblem$modfun)),
                     list(parms=parms,
                          derparms=derparms,
                          code=subproblem$code,
                          bolus=subproblem$bolus,
                          infusion=subproblem$infusion,
                          xdata=subproblem$data$xdata,
                          covdata=subproblem$cov,
                          check=check,
                          options=subproblem$solver.options))
  } else {
    fpred <- do.call(eval(parse(text=subproblem$modfun)),
                     list(parms=parms,
                          derparms=derparms,
                          code=subproblem$code,
                          bolus=subproblem$bolus,
                          infusion=subproblem$infusion,
                          xdata=subproblem$data$xdata,
                          covdata=subproblem$cov,
                          check=check))
  }
  
  # Convert fpred to matrix if is a vector
  if (size(fpred,1)==1)
    fpred <- matrix(fpred,nrow=1)
  
  # Check if fpred contains imaginary numbers
  if (any(is.complex(fpred)))
    cat(' Imaginary number(s) predicted.\n')
  
  # Evaluate the variance function varfun
  weight <- weighting(parms=parms,
                      derparms=derparms,
                      codevar=subproblem$code$var,
                      y=fpred,
                      xdata=subproblem$data$xdata,
                      check=check)
  
  varargout <- list(f=rbind(subproblem$data$xdata,fpred),
                    weight=rbind(subproblem$data$xdata,weight))
  
  return(varargout)

}

