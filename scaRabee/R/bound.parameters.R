
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

bound.parameters <- function(x=NULL,lb=NULL,ub=NULL){

  if (any(dim(x)!=dim(lb)) | any(dim(x)!=dim(ub))){
    stop('  x, lb, and ub arguments do not have the same dimensions.')
  }
  
  newx <- x
  
  for (i in 1:length(x)) {
    if (x[i] < lb[i]){
      newx[i] <- lb[i]
    }
    if (x[i] > ub[i]) {
      newx[i] <- ub[i]
    }
  }
  
  return(newx)
  
}

