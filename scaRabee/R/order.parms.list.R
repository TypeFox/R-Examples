
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

order.parms.list <- function(x=NULL){
  
  index <- list(mod=NULL,
                var=NULL,
                lag=NULL,
                ic =NULL)
  
  for (i in 1:length(x$type)){
    
    if (x$type[i]=='P'){
      index$mod <- c(index$mod,i)
    }
    if (x$type[i]=='V'){
      index$var <- c(index$var,i)
    }
    if (x$type[i]=='L'){
      index$lag <- c(index$lag,i)
    }
    if (x$type[i]=='IC'){
      index$ic <- c(index$ic,i)
    }
  }
  
  ordered <- list(names=c(x$names[index$mod],x$names[index$lag],
                          x$names[index$var],x$names[index$ic]),
                  type =c(x$type[index$mod],x$type[index$lag],
                          x$type[index$var],x$type[index$ic]),
                  value=c(x$value[index$mod],x$value[index$lag],
                          x$value[index$var],x$value[index$ic]),
                  isfix=c(x$isfix[index$mod],x$isfix[index$lag],
                          x$isfix[index$var],x$isfix[index$ic]),
                  lb   =c(x$lb[index$mod],x$lb[index$lag],
                          x$lb[index$var],x$lb[index$ic]),
                  ub   =c(x$ub[index$mod],x$ub[index$lag],
                          x$ub[index$var],x$ub[index$ic]))

  return(as.data.frame(ordered))

}

