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



#' This function creates and trains a self-organizing map (SOM).
#' SOMs are neural networks with one hidden layer. 
#' The network structure is similar to LVQ, but the method is unsupervised 
#' and uses a notion of neighborhood between the units. 
#' The general idea is that the map develops by itself a notion of similarity among 
#' the input and represents this as spatial nearness on the map.
#' Every hidden unit represents a prototype. The goal of learning is to
#' distribute the prototypes in the feature space such that the (probability 
#' density of the) input is represented well.
#' SOMs are usually built with 1d, 2d quadratic, 2d hexagonal, or 3d 
#' neighborhood, so that they can be visualized straightforwardly.
#' The SOM implemented in SNNS has a 2d quadratic neighborhood.
#'  
#' As the computation of this function might be slow if many patterns are involved,
#' much of its output is made switchable (see comments on return values).  
#' 
#' Internally, this function uses the initialization function \code{Kohonen_Weights_v3.2},
#' the learning function \code{Kohonen}, and the update function \code{Kohonen_Order} of
#' SNNS. 
#' 
#' @title Create and train a self-organizing map (SOM)
#' @references 
#' Kohonen, T. (1988), Self-organization and associative memory, Vol. 8, Springer-Verlag.
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
som <- function(x, ...) UseMethod("som")


#' @param x a matrix with training inputs for the network
#' @param mapX the x dimension of the som
#' @param mapY the y dimension of the som
#' @param maxit maximum of iterations to learn
#' @param initFuncParams the parameters for the initialization function
#' @param learnFuncParams the parameters for the learning function
#' @param updateFuncParams the parameters for the update function
#' @param shufflePatterns should the patterns be shuffled?
#' @param calculateMap should the som be calculated?
#' @param calculateActMaps should the activation maps be calculated?
#' @param calculateSpanningTree should the SNNS kernel algorithm for generating a spanning tree be applied?
#' @param saveWinnersPerPattern should a list with the winners for every pattern be saved?
#' @param targets optional target classes of the patterns
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object. Depending on which calculation flags are 
#' switched on, the som generates some special members:
#' \item{map}{the som. For each unit, the amount of patterns where this unit won is given.}
#' \item{componentMaps}{a map for every input component, showing where in the map this component leads to high activation.}
#' \item{actMaps}{a list containing for each pattern its activation map, i.e. all unit activations.
#' The \code{actMaps} are an intermediary result, from which all other results can be computed. This list can be very long,
#' so normally it won't be saved.}
#' \item{winnersPerPattern}{a vector where for each pattern the number of the winning unit is given. Also, an intermediary result
#'  that normally won't be saved.}
#' \item{labeledUnits}{a matrix which -- if the \code{targets} parameter is given -- contains for each unit (rows) and each class 
#' present in the \code{targets} (columns), the amount of patterns of the class where the unit has won. From the \code{labeledUnits}, 
#' the \code{labeledMap} can be computed, e.g. by voting of the class labels for the final label of the unit.}  
#' \item{labeledMap}{a labeled som that is computed from \code{labeledUnits} using \code{\link{decodeClassLabels}}.}
#' \item{spanningTree}{the result of the original SNNS function to calculate the map. For each unit, the last pattern where this unit won is present.
#' As the other results are more informative, the spanning tree is only interesting, if the other functions are too slow or if the original SNNS
#' implementation is needed.}
#' @export
# @S3method som default
#' @method som default
#' @rdname som
#' @examples 
#' \dontrun{demo(som_iris)}
#' \dontrun{demo(som_cubeSnnsR)}
#' 
#' 
#' data(iris)
#' inputs <- normalizeData(iris[,1:4], "norm")
#' 
#' model <- som(inputs, mapX=16, mapY=16, maxit=500,  
#'                 calculateActMaps=TRUE, targets=iris[,5])
#' 
#' par(mfrow=c(3,3))
#' for(i in 1:ncol(inputs)) plotActMap(model$componentMaps[[i]], 
#'                                        col=rev(topo.colors(12)))
#' 
#' plotActMap(model$map, col=rev(heat.colors(12)))
#' plotActMap(log(model$map+1), col=rev(heat.colors(12)))
#' persp(1:model$archParams$mapX, 1:model$archParams$mapY, log(model$map+1), 
#'      theta = 30, phi = 30, expand = 0.5, col = "lightblue")
#' 
#' plotActMap(model$labeledMap)
#' 
#' model$componentMaps
#' model$labeledUnits
#' model$map
#' 
#' names(model)
som.default <- function(x, mapX=16, mapY=16, maxit=100, 
    initFuncParams = c(1.0,  -1.0), 
    learnFuncParams=c(0.5, mapX/2, 0.8, 0.8, mapX), 
    updateFuncParams = c(0.0, 0.0, 1.0),  
    shufflePatterns = TRUE,
    calculateMap=TRUE, calculateActMaps=FALSE, calculateSpanningTree=FALSE, saveWinnersPerPattern=FALSE, targets=NULL, ...) {

  #for the som, no other init/learn/update functions are feasible, so they are not present as
  #parameters..
  initFunc = "Kohonen_Weights_v3.2"
  learnFunc = "Kohonen"
  updateFunc = "Kohonen_Order"
  
  x <- as.matrix(x)

  nInputs <- dim(x)[2L]
  nOutputs <- mapX*mapY
  
  snns <- rsnnsObjectFactory(subclass=c("som"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)
  
  snnsObject <- snns$snnsObject
  
  snnsObject$setLearnFunc(snns$learnFunc)
  snnsObject$setUpdateFunc(snns$updateFunc)
  snnsObject$setInitialisationFunc(snns$initFunc)
  
  snnsObject$setUnitDefaults(0,0,3,0,1,'Act_Logistic','Out_Identity')
  snnsObject$kohonen_createNet(mapX,mapY,nInputs,mapX*mapY)
    
  snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
  snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_Euclid")
  
  patSet <- snnsObject$createPatSet(x)
  
  snnsObject$initializeNet(snns$initFuncParams)
  snnsObject$shufflePatterns(snns$shufflePatterns)
  snnsObject$DefTrainSubPat()
  
#parameters <- c(0.855, 9.0, 0.9, 0.9, 16.0)
#parameters <- c(0.855, mapX/2, 0.9, 0.9, mapX)
  #parameters <- c(0.5, mapX/2, 0.8, 0.8, mapX)
  #maxit <- 1000
    
  for(i in 1:maxit) {
    res <- snnsObject$learnAllPatterns(snns$learnFuncParams)
  }
  
  snns$archParams <- list(mapX=mapX, mapY=mapY)
  
  snns$nOutputs <- nOutputs
  
  if(calculateMap) {

    mapVec <- snnsObject$somPredictCurrPatSetWinners(updateFuncParams, 
                                  saveWinnersPerPattern=saveWinnersPerPattern, targets=targets)
    snns$map <- vectorToActMap(mapVec$nWinnersPerUnit, nrow=mapX)
    snns$winnersPerPattern <- mapVec$winnersPerPattern
    
    if(!is.null(targets)) {
      
      snns$labeledUnits <- mapVec$labeledUnits
      labeledMap <- encodeClassLabels(mapVec$labeledUnits, method="WTA", l=0,h=0)
      snns$labeledMap <- vectorToActMap(labeledMap, mapX)      
        
    } else {
      snns$labeledUnits <- NULL
      snns$labeledMap <- NULL          
    }
    
#    if(!is.null(targets)) {
#      labels <- as.numeric(targets)
#      labeledSpanningTree <- spanningTree
#      for(i in 1:nrow(spanningTree))
#        for(j in 1:ncol(spanningTree))
#        {
#          if(spanningTree[i,j] != 0)  {
#            labeledSpanningTree[i,j] <- labels[spanningTree[i,j]]      
#          } 
#        }
#      snns$labeledSpanningTree <- labeledSpanningTree
#    } else {
#      snns$labeledSpanningTree <- NULL
#    }
    
  } else {
    snns$map <- NULL
    snns$winnersPerPattern <- NULL
    snns$labeledUnits <- NULL
    snns$labeledMap <- NULL
  }

  if(calculateActMaps) {
    
    actMat <- snnsObject$predictCurrPatSet("som", updateFuncParams)
    snns$actMaps <- matrixToActMapList(actMat, nrow=mapX)
    
  } else {
    snns$actMaps <- NULL
  }
  
  if(calculateSpanningTree) {

    spanningTree <- snnsObject$somPredictCurrPatSetWinnersSpanTree()
    snns$spanningTree <- vectorToActMap(spanningTree, nrow=mapX)
    
  } else {
    snns$spanningTree <- NULL
  }

  #calculate component maps
  compMat <- snnsObject$somPredictComponentMaps(updateFuncParams)  
  snns$componentMaps <- matrixToActMapList(compMat, nrow=mapX)
  
  snnsObject$deletePatSet(patSet$set_no)
  
  snns
}

