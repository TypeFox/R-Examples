
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

input.scaling <- function(parms=NULL,
                          derparms=NULL,
                          codescale=NULL,
                          ic=NULL,
                          check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(codescale))
      stop('codescale argument is NULL.')
    
    if (!is.character(codescale))
      stop('codescale argument must be an object of class \'character\'')
    
    if (length(codescale)>1)
      stop('codescale argument must have a length of 1.')
  }
  
  # Evaluation of codescale
  if (is.null(derparms)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,derparms)
  }
  
  scale <- with(cparms,{
    
    eval(parse(text=codescale))
    
    return(scale)
  })
  
  # Check scale object
  if (check){
    if (is.null(scale))
      stop('scale object created in the $SCALE block cannot be NULL.')
    
    if (size(scale,1)!=1)
      stop('scale must be a scalar or a vector.')
  }
  
  # Expand scale if it is a scalar
  if (size(scale,2)==1){
    scale <- rep(scale,size(ic,2))
  } else {
    if (check){
      if (size(scale,2)!=size(ic,2))
        stop(paste('scale must be a scalar or have the same',
                   'dimension as the system of\n  differential equations.'))
    }
  }
  
  return(scale)
  
}

