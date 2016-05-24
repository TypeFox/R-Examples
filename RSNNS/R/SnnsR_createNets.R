#############################################################################
#
#   This file is part of the R package "RSNNS".
#
#   Author: Christoph Bergmeir
#   Supervisor: José M. Benítez
#   Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
# 
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
# 
#   You should have received a copy of the GNU Library General Public License
#   along with this library; see the file COPYING.LIB.  If not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#   Boston, MA 02110-1301, USA.
#
#############################################################################

# This is an SnnsR low-level function. You may want to have a look 
# at \code{\link{SnnsR-class}} to find out how to properly use it.
# @seealso \code{\link{SnnsR-class}}

# @include SnnsWrapperFunctions.R
#NULL

#' This function creates a layered network in the given SnnsR object.
#' This is an SnnsR low-level function. You may want to have a look 
#' at \code{\link{SnnsR-class}} to find out how to properly use it.
#'  
#' @title Create a layered network
#' @param unitsPerLayer a vector of integers that represents the number of units in each layer, including input and output layer
#' @param fullyConnectedFeedForward if \code{TRUE}, the network is fully connected as a feed-forward network. If \code{FALSE}, 
#' no connections are created
#' @param iNames names of input units
#' @param oNames names of output units
#' @rdname SnnsRObject-createNet
#' @name SnnsRObject$createNet
#' @usage \S4method{createNet}{SnnsR}(unitsPerLayer, fullyConnectedFeedForward = TRUE, iNames = NULL, oNames = NULL)
#' @aliases createNet,SnnsR-method SnnsR__createNet
#' @seealso \code{\link{SnnsR-class}}
#' @examples
#' obj1 <- SnnsRObjectFactory()
#' obj1$createNet(c(2,2), FALSE)
#' obj1$getUnitDefinitions()
#' 
#' obj2 <- SnnsRObjectFactory()
#' obj2$createNet(c(8,5,5,2), TRUE)
#' obj2$getUnitDefinitions()
SnnsR__createNet <- function(snnsObject, unitsPerLayer, fullyConnectedFeedForward = TRUE, iNames = NULL, oNames = NULL) {
  
  if(length(unitsPerLayer) < 2) stop("At least 2 layers have to be specified")
  
  layers <- list()
  currLayer <- 1
  
  nInputs <- unitsPerLayer[currLayer]
  layers[[currLayer]] <- vector()

  if( is.null(iNames) )
    iNames <- paste("Input", 1:nInputs, sep="_")
  else {
    iNames <- paste("Input", as.character(iNames), sep="_")
    stopifnot( length(iNames) == nInputs )
  }
  
  for(i in 1:nInputs) {
    
    num <- snnsObject$createDefaultUnit()
    layers[[currLayer]][i] <- num
    
    snnsObject$setUnitName(num, iNames[[i]])
    
    snnsObject$setUnitTType(num, resolveSnnsRDefine("topologicalUnitTypes","UNIT_INPUT"))
    
    snnsObject$setUnitPosition(num, i, 0, 0)
    
  }
  
  currLayer <- currLayer + 1
  
  
  for(k in seq(length=(length(unitsPerLayer)-2)))  {
    
    nHidden <- unitsPerLayer[currLayer]
    layers[[currLayer]] <- vector()
    
    for(i in 1:nHidden) {
      
      num <- snnsObject$createDefaultUnit()
      layers[[currLayer]][i] <- num
      
      snnsObject$setUnitName(num,paste("Hidden_",currLayer,"_",i,sep=""))
     
      snnsObject$setUnitTType(num, resolveSnnsRDefine("topologicalUnitTypes","UNIT_HIDDEN"))
      
      snnsObject$setUnitPosition(num, i, (currLayer-1)*2, 0)
      
      if(fullyConnectedFeedForward)  {
        
        snnsObject$setCurrentUnit(num)
        
        for(j in layers[[currLayer-1]])  {
          snnsObject$createLink(j,0);
        }   
      }
    }
    
    currLayer <- currLayer + 1
  }
  
  nOutputs <- unitsPerLayer[currLayer]
  layers[[currLayer]] <- vector()

  
  if( is.null(oNames) )
    oNames <- paste("Output", 1:nOutputs, sep="_")
  else {
    oNames <- paste("Output", as.character(oNames), sep="_")
    stopifnot( length(oNames) == nOutputs )
  }
  
  for(i in 1:nOutputs) {
    
    num <- snnsObject$createDefaultUnit()
    layers[[currLayer]][i] <- num
    
    snnsObject$setUnitName(num, oNames[[i]])
    
    snnsObject$setUnitTType(num, resolveSnnsRDefine("topologicalUnitTypes","UNIT_OUTPUT"))
    
    snnsObject$setUnitPosition(num, i, (currLayer-1)*2, 0)
    
    if(fullyConnectedFeedForward)  {

      snnsObject$setCurrentUnit(num)
      
      for(j in layers[[currLayer-1]])  {
        snnsObject$createLink(j,0);
      }
    }
  }
} 

