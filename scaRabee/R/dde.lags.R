
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

dde.lags <- function(parms=NULL,
                     derparms=NULL,
                     codelags=NULL,
                     check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(parms))
      stop('parms argument (primary parameters) is NULL.')
    
    if (is.null(codelags))
      stop('codelags argument (derived parameters) is NULL.')
    
    if (!is.character(codelags))
      stop('codelags argument must be an object of class \'character\'')
    
    if (length(codelags)>1)
      stop('codelags argument must have a length of 1.')
  }
  
  # Evaluation of codeparms
  if (is.null(derparms)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,derparms)
  }
  
  lags <- with(cparms,{
    # Get all objects available before evaluation
    objects.start <- ls(all=TRUE)
    
    # Evaluation
    eval(parse(text=codelags))
    
    # Get all objects created by evaluation
    objects.end <- ls(all=TRUE)
    objects.end <- objects.end[-which(objects.end=='objects.start')]
    
    objects.new <- lapply(objects.end[!objects.end%in%objects.start],
                          function(x) eval(parse(text=x)))
    
    names(objects.new) <- objects.end[!objects.end%in%objects.start]
    return(objects.new)
  })
  
  # Check lags object
  if (!is.null(lags))
    lags <- unlist(lags)
  
  # Check for overwritten L parameters
  lparms <- parms[which(attr(parms,'type')=='L')]
  if (length(lparms)>0){
    overwritten <- numeric(0)
    
    if (!is.null(lags))
      overwritten <- which(names(lags)%in%names(lparms))
    
    if (length(overwritten)>0){
      lags <- lags[-overwritten]
    } 
    
    lags <- c(lags, lparms)
  }
  
  if (check && is.null(lags))
    stop('no system delay was provided.')
  
  return(lags)
  
}