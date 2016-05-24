
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

weighting <- function(parms=NULL,
                      derparms=NULL,
                      codevar=NULL,
                      y=NULL,
                      xdata=NULL,
                      check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(codevar))
      stop('codevar argument is NULL.')
    
    if (!is.character(codevar))
      stop('codevar argument must be an object of class \'character\'')
    
    if (length(codevar)>1)
      stop('codevar argument must have a length of 1.')
  }
  
  # Evaluation of codevar
  if (is.null(derparms)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,derparms)
  }
  
  # Get dimension of xdata
  ntime <- length(xdata)
  
  v <- with(cparms,{
      
      eval(parse(text=codevar))
      
      return(v)
      
  })
   
  # Check dimension of v
  if (check){
    if (sum(size(v)==size(y))!=2)
      stop(paste('the dimensions of the variance matrix should not differ from',
                 'the dimensions of\n  the predictions matrix. Please, check',
                 'your code.'))
  }
  
  return(v)
  
}