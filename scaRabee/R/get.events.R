
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

get.events <- function(bolus=NULL, scale=NULL){
  
  if (is.null(bolus) || dim(bolus)[1]==0)
    return(data.frame(var=character(0),
                      time=numeric(0),
                      value=numeric(0),
                      method=character(0)))
  
  events <- data.frame(var=paste('a',bolus$CMT,sep=''),
                       time=bolus$TIME,
                       value=bolus$AMT,
                       method=rep('add',dim(bolus)[1]))
                     
  
  if (!is.null(scale)){
    events$value <- events$value/scale[bolus$CMT]
  }
  
  return(events)
  
}