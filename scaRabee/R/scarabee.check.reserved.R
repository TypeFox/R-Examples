
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

scarabee.check.reserved <- function(names=NULL, covnames=NULL){
  
  reserved <- c(paste('y',1:25,sep=''),
                'x', 'dosing', 'xdata', 'covdata', 'issim', 'parms', 'tspan', 
                't', 'time', 'f', 'scale', 'dydt', 'delays')
  
  if (!is.null(names)){
    if (any(reserved%in%names)){
      hits <- paste(reserved[reserved%in%names], collapse=', ')
      stop(paste('model parameters cannot use reserved names.',
                 'Please, rename the following parameter(s):\n  ',hits))
    }
  }
  
  if (!(is.null(names) | is.null(covnames))){
    if (any(toupper(covnames)%in%names)){
      hits <- paste(toupper(covnames)[toupper(covnames)%in%names], collapse=', ')
      stop(paste('model parameters cannot use covariate names.',
                 'Please, rename the following parameter(s):\n  ',hits))
    }
  }
}
