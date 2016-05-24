
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

compute.secondary <- function(subproblem=NULL,x=NULL){
  
  secondary <- function(parms=NULL,codesec=NULL){
    
    # Input validation
    if (is.null(parms))
      stop('parms argument (primary parameters) is NULL.')
    
    if (is.null(codesec))
      stop('codesec argument (secondary parameters) is NULL.')
    
    if (!is.character(codesec))
      stop('codesec argument must be an object of class \'character\'')
    
    if (length(codesec)>1)
      stop('codesec argument must have a length of 1.')
    
    # Evaluation of codeparms
    objects.new <- with(as.list(parms),{
      # Get all objects available before evaluation
      objects.start <- ls(all=TRUE)
      
      # Evaluation
      eval(parse(text=codesec))
      
      # Get all objects created by evaluation
      objects.end <- ls(all=TRUE)
      objects.end <- objects.end[-which(objects.end=='objects.start')]
      
      objects.new <- lapply(objects.end[!objects.end%in%objects.start],
                            function(x) eval(parse(text=x)))
      
      names(objects.new) <- objects.end[!objects.end%in%objects.start]
      return(objects.new)
    })
    
    secparms <- unlist(objects.new)
    
    return(secparms)
    
  }
  
  # Copies subproblem.init in newparam and replaces the new estimates in newparam
  newparam <- initparam <- subproblem$init
  npar     <- length(newparam$names)
  estindex <- 1
  for (i in 1:npar){
    if (newparam$isfix[i]==0){
      newparam$value[i] <- x[estindex]
      estindex <- estindex + 1
    }
  }
  
  # Retrieve initial primary parameters
  parms <- c(get.parms.data(x=initparam,which='value',type='P'),
             get.parms.data(x=initparam,which='value',type='L'),
             get.parms.data(x=initparam,which='value',type='IC'),
             get.parms.data(x=initparam,which='value',type='V'))
  names(parms) <- c(get.parms.data(x=initparam,which='names',type='P'),
                    get.parms.data(x=initparam,which='names',type='L'),
                    get.parms.data(x=initparam,which='names',type='IC'),
                    get.parms.data(x=initparam,which='names',type='V'))
  
  if (!is.null(subproblem$code$sec)){
    # Computes initial value of secondary parameters
    init <- secondary(parms=parms,
                      codesec=subproblem$code$sec)
    
    # Retrieve final primary parameters
    parms <- c(get.parms.data(x=newparam,which='value',type='P'),
               get.parms.data(x=newparam,which='value',type='L'),
               get.parms.data(x=newparam,which='value',type='IC'),
               get.parms.data(x=newparam,which='value',type='V'))
    names(parms) <- c(get.parms.data(x=newparam,which='names',type='P'),
                      get.parms.data(x=newparam,which='names',type='L'),
                      get.parms.data(x=newparam,which='names',type='IC'),
                      get.parms.data(x=newparam,which='names',type='V'))
    
    # Computes secondary parameter estimates
    estimates <- secondary(parms=parms,
                           codesec=subproblem$code$sec)
    
    varargout <- list(init=init,estimates=estimates,names=names(init))
  } else {
    varargout <- list(init=NULL,estimates=NULL,names=NULL)
  }
  
  return(varargout)
  
}

