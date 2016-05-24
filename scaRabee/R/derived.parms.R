
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

derived.parms <- function(parms=NULL,
                          covdata=NULL,
                          codederiv=NULL,
                          check=FALSE){
  
  # Input validation
  if (check){
    if (is.null(parms))
      stop('parms argument (primary parameters) is NULL.')
    
    if (is.null(codederiv))
      stop('codederiv argument (derived parameters) is NULL.')
    
    if (!is.character(codederiv))
      stop('codederiv argument must be an object of class \'character\'')
    
    if (length(codederiv)>1)
      stop('codederiv argument must have a length of 1.')
  }
  
  # Prevent modification of parameters and covariates
  redef.parms <- any(unlist(sapply(
    names(parms),
    function(x,codederiv) {grep(paste(x,"[[:space:]]*[=<]",sep=""),codederiv)},
    codederiv)))
  
  redef.covdata <- any(unlist(sapply(
    names(covdata),
    function(x,codederiv) {grep(paste(x,"[[:space:]]*[=<]",sep=""),codederiv)},
    codederiv)))
  
  if (redef.parms)
    stop('Parameters cannot be modified in $DERIVED.')
  
  if (redef.covdata)
    stop('Covariates cannot be modified in $DERIVED.')
  
  # Evaluation of codederiv
  if (is.null(covdata)){
    cparms <- as.list(parms)
  } else {
    cparms <- c(parms,covdata)
  }
  
  objects.new <- with(cparms,{
    # Get all objects available before evaluation
    objects.start <- ls(all=TRUE)
    
    # Evaluation
    eval(parse(text=codederiv))
    
    # Get all objects created by evaluation
    objects.end <- ls(all=TRUE)
    objects.end <- objects.end[-which(objects.end=='objects.start')]
    
    objects.new <- lapply(objects.end[!objects.end%in%objects.start],
      function(x) eval(parse(text=x)))
    
    names(objects.new) <- objects.end[!objects.end%in%objects.start]
    return(objects.new)
  })
  
  return(objects.new)
  
}
