
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

init.update <- function(a=NULL,
                        t=NULL,
                        dosing=NULL,
                        scale=NULL){
  
  # Set init to a
  if (names(a)[1]=='time'){
    init <- a[-1]
  } else {
    init <- a
  }
  
  # Subset dosing for event occuring at time t
  bolus <- dosing[dosing[,1]==t,,drop=FALSE]
  
  # Update init
  if (any(bolus[,3]>0)){
    for (i in 1:size(bolus,1)) {
      init[bolus[i,2]] <- init[bolus[i,2]] + bolus[i,3]/scale[bolus[i,2]]
    }
  }
  
  return(init)
  
}
