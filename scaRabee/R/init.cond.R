
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

init.cond <- function(parms=NULL,
                      derparms=NULL,
                      codeic=NULL,
                      dosing=NULL,
                      check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(dosing))
      stop('dosing argument is NULL.')
    
    if (is.null(codeic))
      stop('codeic argument is NULL.')
    
    if (!is.character(codeic))
      stop('codeic argument must be an object of class \'character\'')
    
    if (length(codeic)>1)
      stop('codeic argument must have a length of 1.')
  }
  
  # Evaluation of codeic
  if (is.null(derparms)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,derparms)
  }
  
  init <- with(cparms,{
    
    eval(parse(text=codeic))
    names(init) <- paste('a',1:length(init),sep='')
    return(init)
    
  })
  
  # Check init object
  if (check){
    if (is.null(init))
      stop('init object created in the $IC block cannot be NULL.')
    
    nstate <- size(init,2)
    
    if (any(is.na(match(dosing[,2],c(1:nstate)))))
      stop(paste('one or more inputs are assigned to a state that is not ',
                 'defined in the system\n  of differential equations.',sep=''))
  }
  
  return(init)
  
}