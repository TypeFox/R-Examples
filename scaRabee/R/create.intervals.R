
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

create.intervals <- function(xdata=NULL,bolus=NULL,infusion=NULL){
  
  # Trim bolus of any duplicates time records
  bolus <- bolus[match(unique(bolus[,1]),bolus[,1]),]
  
  # Determines the min time to start the integration. The maximum has to be
  # the maximum observation time.
  if (size(bolus,1)>0) {
    mintime <- min(c(transpose(bolus[,1]),xdata))
  } else {
    if (size(infusion,1)>0) {
      mintime <- min(c(transpose(infusion[,1]),xdata))
    } else {
      mintime <- min(xdata)
    }
  }
  maxtime <- max(xdata)
  
  # Trims bolus from doses over maxtime
  bolus <- bolus[bolus[,1]<=maxtime,]
  
  # Initializes intervals
  if (size(bolus,1)>0) {
    if (bolus[1,1]==mintime){   # integration starts at first dose
      intervals <- c()
    } else {                    # integration starts before first dose
      intervals <- rbind(mintime,bolus[1,1])
    }
      
    # Buils intervals
    ndose <- size(bolus,1)
    for (i in 1:ndose){
      if (i!=ndose){                 # intermediate dose
        tmp <- rbind(bolus[i,1],bolus[i+1,1])
      } else {                       # last dose
        if (bolus[i,1]==maxtime){    # time of last dose = time of last observation
          tmp <- NULL
        } else {                     # otherwise
          tmp <- rbind(bolus[i,1],maxtime)
        }
      }
      intervals <- cbind(intervals,tmp)
    }
  } else {
    intervals <- rbind(mintime,maxtime)
  }
  row.names(intervals) <- c('start.time','end.time')
  
  return(intervals)
  
}

