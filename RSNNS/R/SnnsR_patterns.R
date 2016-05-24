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


# @include SnnsRObjectFactory.R
#NULL

#' SnnsR low-level function to create a pattern set in the SNNS kernel from 
#' the values given, so that they are available in the SNNS kernel for use.
#' 
#' @title Create a pattern set
#' @param inputs the input values
#' @param targets the target values
#' @return a list with elements \code{err} and \code{set_no}. The latter one identifies the pattern set within the \code{\link{SnnsR-class}} object 
#' @rdname SnnsRObject-createPatSet
#' @name SnnsRObject$createPatSet
#' @usage \S4method{createPatSet}{SnnsR}(inputs, targets)
#' @aliases createPatSet,SnnsR-method SnnsR__createPatSet
SnnsR__createPatSet <- function(snnsObject, inputs, targets) {

  #sort is necessary to fix problems with pruned networks, where the order may change
  iUnits <- sort(snnsObject$getAllInputUnits())
  oUnits <- sort(snnsObject$getAllOutputUnits())
  
  x <- as.matrix(inputs)
  nInputs <- ncol(x)
  if (nInputs != length(iUnits)) 
    stop(paste("number of input data columns ", nInputs ," does not match number of input neurons ", length(iUnits) ,sep=""))
  
  if(!missing(targets)){
    y <- as.matrix(targets)
    nOutputs <- ncol(y)
    if (nOutputs != length(oUnits)) 
      stop(paste("number of output data columns ", nOutputs ," does not match number of output neurons ", length(oUnits) ,sep=""))
    
    if(nrow(x) != nrow(y)) 
      stop(paste("input value rows ",nrow(x)," not same as output value rows ",nrow(y),sep=""))
  }
  
  patSet <- snnsObject$allocNewPatternSet()


  
  for(i in 1:nrow(x)) {
    for(j in 1:nInputs)  {
      #snnsObject$setUnitActivation(iUnits[(nInputs+1)-j], x[i,j])
      snnsObject$setUnitActivation(iUnits[j], x[i,j])
    }
    
    if(!missing(targets) && length(targets) != 0) {  
      for(j in 1:nOutputs)  {
        #snnsObject$setUnitActivation(oUnits[(nOutputs+1)-j], y[i,j])
        snnsObject$setUnitActivation(oUnits[j], y[i,j])
      }
    }
    snnsObject$newPattern()
  }

  snnsObject$setCurrPatSet(patSet$set_no)
  
  return(patSet)
}


#' SnnsR low-level function for generic prediction with a trained net.
#' 
#' @title Predict values with a trained net
#' @param units the units that define the output
#' @param updateFuncParams the parameters for the update function (the function has to be already set)
#' @return the predicted values
#' @rdname SnnsRObject-genericPredictCurrPatSet
#' @name SnnsRObject$genericPredictCurrPatSet
#' @usage \S4method{genericPredictCurrPatSet}{SnnsR}(units, updateFuncParams=c(0.0))
#' @aliases genericPredictCurrPatSet,SnnsR-method SnnsR__genericPredictCurrPatSet
SnnsR__genericPredictCurrPatSet <- function(snnsObject, units, updateFuncParams=c(0.0))  {
  
  noOfPatterns <- snnsObject$getNoOfPatterns()
  
  predictions <- matrix(nrow= noOfPatterns, ncol=length(units))
  
  snnsObject$DefTrainSubPat()  
  
  for(currentPattern in 1:noOfPatterns)  {
    
    snnsObject$setPatternNo(currentPattern)
    
    snnsObject$showPattern(resolveSnnsRDefine("patternUpdateModes","OUTPUT_NOTHING"))
    
    snnsObject$updateNet(updateFuncParams)
    
    for(i in 1:length(units)) {
      predictions[currentPattern,i] <- snnsObject$getUnitOutput(units[i])
    }
    
  }
  
  return(predictions)
} 


#' SnnsR low-level function to get a list of output units of a net.
#' 
#' Depending on the network architecture, output is present in hidden units, in
#' output units, etc. In some network types, the output units have a certain
#' name prefix in SNNS. This function finds the output units according to
#' certain network types. The type is specified by \code{outputMethod}. If the
#' given \code{outputMethod} is unknown, the function defaults to "output".
#' 
#' @title Get a list of output units of a net
#' @param outputMethod a string defining the output method of the net. Possible values are: "art1", "art2", "artmap", "assoz", "som", "output".
#' @return a list of numbers identifying the units
#' @rdname SnnsRObject-whereAreResults
#' @name SnnsRObject$whereAreResults
#' @usage \S4method{whereAreResults}{SnnsR}(outputMethod="output")
#' @aliases whereAreResults,SnnsR-method SnnsR__whereAreResults
SnnsR__whereAreResults <- function(snnsObject, outputMethod="output") {
  
  units <- NULL
  #outputMethod <- "art1"
  
  if(outputMethod == "art1") {  
    # in the ART1 network, the units that represent the output patterns are named cmp1, cmp2, ...
    #units <- snnsObject$getUnitsByName("cmp")
    units <- snnsObject$getUnitsByName("rec")
    
  } else if(outputMethod == "art2") {
    
    #unitsX <- snnsObject$getUnitsByName("x")
    #unitsQ <- snnsObject$getUnitsByName("q")
    #units <- c(unitsX, unitsQ) 
    
    units <- snnsObject$getUnitsByName("rec")
    
  } else if(outputMethod == "artmap") {
    
    units <- snnsObject$getUnitsByName("map")
        
  } else if(outputMethod=="assoz") {
    
    units <- snnsObject$getAllHiddenUnits()
    
  } else if(outputMethod=="som") {
    
    units <- snnsObject$getAllHiddenUnits()
    
  } else { #if(outputMethod=="reg_class") {
    
    units <- snnsObject$getAllOutputUnits()
    
  } #else if(outputMethod=="assoz") {  }
  
  units
}

#' SnnsR low-level function to predict values with a trained net.
#' 
#' This function has to be used embedded in a step of loading and afterwards 
#' removing the patterns into the \code{\link{SnnsR-class}} object. As SNNS only supports 2 pattern sets
#' in parallel, removing unneeded pattern sets is quite important.
#' 
#' @title Predict values with a trained net
#' @param outputMethod is passed to \link{SnnsRObject$whereAreResults}
#' @param updateFuncParams parameters passed to the networks update function
#' @return the predicted values
#' @rdname SnnsRObject-predictCurrPatSet
#' @name SnnsRObject$predictCurrPatSet
#' @usage \S4method{predictCurrPatSet}{SnnsR}(outputMethod="reg_class", updateFuncParams=c(0.0))
#' @aliases predictCurrPatSet,SnnsR-method SnnsR__predictCurrPatSet
SnnsR__predictCurrPatSet <- function(snnsObject, outputMethod="reg_class", updateFuncParams=c(0.0))  {
  
  units <- snnsObject$whereAreResults(outputMethod)
  predictions <- snnsObject$genericPredictCurrPatSet(units, updateFuncParams)
  predictions
}


#' SnnsR low-level function to calculate the som component maps.
#' 
#' @title Calculate the som component maps
#' @param updateFuncParams parameters passed to the networks update function
#' @return a matrix containing all componant maps as 1d vectors
#' @rdname SnnsRObject-somPredictComponentMaps
#' @name SnnsRObject$somPredictComponentMaps
#' @usage \S4method{somPredictComponentMaps}{SnnsR}(updateFuncParams=c(0.0, 0.0, 1.0))
#' @aliases somPredictComponentMaps,SnnsR-method SnnsR__somPredictComponentMaps
#' @seealso \code{\link{som}}
SnnsR__somPredictComponentMaps <- function(snnsObject, updateFuncParams=c(0.0, 0.0, 1.0))  {
  
  snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_Component")
  
  nInputs <- snnsObject$getNoOfInputUnits()
  units <- snnsObject$getAllHiddenUnits()
  
  predictions <- matrix(nrow= nInputs, ncol=length(units))
  
  
  for(input in 1:nInputs)  {
  
    #parameter has to be reversed to get same order as in SNNS gui.. TODO: why?
    snnsObject$kohonen_SetExtraParameter((nInputs+1) - input)
    snnsObject$updateNet(updateFuncParams)
    
    for(i in 1:length(units)) {
      predictions[input,i] <- snnsObject$getUnitOutput(units[i])
      #predictions[input,i] <- snnsObject$getUnitValueA(units[i])
    }
    
  }
  
  snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_Euclid")
  
  return(predictions)
} 


#' SnnsR low-level function to get most of the relevant results from a SOM.
#'  
#' @title Get most of the relevant results from a som
#' @param updateFuncParams parameters passed to the networks update function
#' @param saveWinnersPerPattern should a list with the winners for every pattern be saved?
#' @param targets optional target classes of the patterns
#' @return a list with three elements:
#' \item{nWinnersPerUnit}{For each unit, the amount of patterns where this unit won is given. So, this is a 1d vector representing the normal version of the som.}
#' \item{winnersPerPattern}{a vector where for each pattern the number of the winning unit is given. This is an intermediary result
#'  that normally won't be saved.}
#' \item{labeledUnits}{a matrix which -- if the \code{targets} parameter is given -- contains for each unit (rows) and each class 
#' present in the \code{targets} (columns), the amount of patterns of the class where the unit has won. From the \code{labeledUnits}, 
#' the \code{labeledMap} can be computed, e.g. by voting of the class labels for the final label of the unit.}
#' @rdname SnnsRObject-somPredictCurrPatSetWinners
#' @name SnnsRObject$somPredictCurrPatSetWinners
#' @usage \S4method{somPredictCurrPatSetWinners}{SnnsR}(updateFuncParams=c(0.0, 0.0, 1.0), 
#' saveWinnersPerPattern=TRUE, targets=NULL)
#' @aliases somPredictCurrPatSetWinners,SnnsR-method SnnsR__somPredictCurrPatSetWinners
#' @seealso \code{\link{som}}
SnnsR__somPredictCurrPatSetWinners <- function(snnsObject, updateFuncParams=c(0.0, 0.0, 1.0), saveWinnersPerPattern=TRUE, targets=NULL)  {
  
  units <- snnsObject$getAllHiddenUnits()
  noOfPatterns <- snnsObject$getNoOfPatterns()
  winners <- snnsObject$somPredictCurrPatSetWinnersC(units, noOfPatterns, updateFuncParams)
  
  map <- seq(0, 0, length=length(units))
  
  for(i in 1:length(winners)) {
    map[winners[i]] <- map[winners[i]] + 1 
  }

  if(!is.null(targets)) {

    classes <- unique(targets)
    numClasses <- 1:length(classes)
    names(numClasses) <- classes
    labeledUnits <- matrix(0, nrow=length(units), ncol=length(classes))
    colnames(labeledUnits) <- classes
    
    for(i in 1:length(winners)) {
      currUnit <- winners[i]
      labeledUnits[currUnit, numClasses[targets[i]]] <- labeledUnits[currUnit, numClasses[targets[i]]] + 1
    }
  } else {
    labeledUnits <- NULL
  }
  
  if(!saveWinnersPerPattern) winners <- NULL
  
  return(list(nWinnersPerUnit=map, winnersPerPattern=winners, labeledUnits=labeledUnits))
 
#  Function is now reimplemented in C++ to be faster..
  
#  units <- snnsObject$getAllHiddenUnits()
#  noOfPatterns <- snnsObject$getNoOfPatterns()
#  
#  winners <- vector()
#  
#  for(currentPattern in 1:noOfPatterns)  {
#   
#    predictions <- vector()
#    
#    snnsObject$setPatternNo(currentPattern)
#    snnsObject$showPattern(resolveSnnsRDefine("patternUpdateModes","OUTPUT_NOTHING"))
#    snnsObject$updateNet(updateFuncParams)
#    
#    for(i in 1:length(units)) {
#      predictions[i] <- snnsObject$getUnitOutput(units[i])
#    }
#    
#    winners <- c(winners, which(predictions == min(predictions, na.rm=TRUE)))
#  }
#
#  map <- seq(0, 0, length=length(units))
#
#  for(i in 1:length(winners)) {
#    map[winners[i]] <- map[winners[i]] + 1 
#  }
#
#  map
  
} 


#' SnnsR low-level function to get the spanning tree of the SOM, This function 
#' calls directly the corresponding SNNS kernel function (the only one available for SOM).
#' Advantage are faster computation, disadvantage is somewhat limited information in
#' the output. 
#' 
#' @title Get the spanning tree of the SOM
#' @return the spanning tree, which is the som, showing for each unit a number identifying 
#' the last pattern for which this unit won. (We note that, also if there are more than 
#' one patterns, only the last one is saved)  
#' @rdname SnnsRObject-somPredictCurrPatSetWinnersSpanTree
#' @name SnnsRObject$somPredictCurrPatSetWinnersSpanTree
#' @usage \S4method{somPredictCurrPatSetWinnersSpanTree}{SnnsR}()
#' @aliases somPredictCurrPatSetWinnersSpanTree,SnnsR-method SnnsR__somPredictCurrPatSetWinnersSpanTree
#' @seealso \code{\link{som}}
SnnsR__somPredictCurrPatSetWinnersSpanTree <- function(snnsObject)  {
  
  units <- snnsObject$getAllHiddenUnits()
  
  noOfPatterns <- snnsObject$getNoOfPatterns()
  
  snnsObject$spanning_tree()
      
  predictions <- vector() 
    
    for(i in 1:length(units)) {
      predictions[i] <- snnsObject$getUnitValueA(units[i])
    }
  
  return(predictions)
} 

