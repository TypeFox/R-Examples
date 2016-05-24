
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

de.output <- function(f=NULL,
                      parms=NULL,
                      derparms=NULL,
                      codeoutput=NULL,
                      dosing=NULL,
                      xdata=NULL,
                      check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(codeoutput))
      stop('codeoutput argument is NULL.')
    
    if (!is.character(codeoutput))
      stop('codeoutput argument must be an object of class \'character\'')
    
    if (length(codeoutput)>1)
      stop('codeoutput argument must have a length of 1.')
  }
  
  # Create time variable from f and drop first row
  times <- f[1,]
  f <- f[-1,,drop=FALSE]
  
  # Evaluation of codeoutput
  if (is.null(derparms)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,derparms)
  }
  
  y <- with(cparms,{
    
    eval(parse(text=codeoutput))
    
    return(y)
  })
    
  if (size(y,1)==1 & size(y,2)>1)
    y <- matrix(y,nrow=1)
  
  return(y)
  
}