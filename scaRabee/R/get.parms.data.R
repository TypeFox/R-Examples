
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

get.parms.data <- function(x=NULL, which=NULL, type=NULL){

  if (!any(which==c('names','type','value','isfix')))
    stop(sprintf('\n  the data frame x does not contain any %s variable.',
                 which))

  if (!any(type==c('P','L','V','IC')))
    stop('the parameter category should be P, L, V, or IC')

  if (size(x,1)==0)
    return(NULL)
  
  mystr <- paste('x$', which,sep='')
  xpar  <- eval(parse(text=mystr))
  xtype <- x$type
  indices <- c()
  for (i in 1:length(xpar)){
    if (xtype[i]==type){
      indices <- c(indices,i)
    }
  }
  param <- xpar[indices]

  return(param)
  
}

