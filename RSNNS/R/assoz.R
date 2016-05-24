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


#' The autoassociative memory performs clustering by finding a prototype to the given input. 
#' The implementation assumes two-dimensional input and output (cf. \code{\link{art1}}).
#' 
#' The default initialization and update functions are the only ones suitable for this kind of 
#' network. The update function takes one parameter, which is the number of iterations that will 
#' be performed. The default of 50 usually does not have to be modified. For learning, \code{RM_delta} 
#' and \code{Hebbian} functions can be used, though the first one usually performs better.
#' 
#' A more detailed description of the theory and the parameters is available from 
#' the SNNS documentation and the other referenced literature. 
#' 
#' @title Create and train an (auto-)associative memory
#' @references 
#' Palm, G. (1980), 'On associative memory', Biological Cybernetics 36, 19-31.
#' 
#' Rojas, R. (1996), Neural networks :a systematic introduction, Springer-Verlag, Berlin.
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' @export
assoz <- function(x, ...) UseMethod("assoz")


#' @param x a matrix with training inputs for the network
#' @param dimX x dimension of inputs and outputs
#' @param dimY y dimension of inputs and outputs
#' @param maxit maximum of iterations to learn
#' @param initFunc the initialization function to use
#' @param initFuncParams the parameters for the initialization function
#' @param learnFunc the learning function to use
#' @param learnFuncParams the parameters for the learning function
#' @param updateFunc the update function to use
#' @param updateFuncParams the parameters for the update function
#' @param shufflePatterns should the patterns be shuffled?
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object. The \code{fitted.values} member contains the 
#' activation patterns for all inputs.
#' @export
# @S3method assoz default
#' @method assoz default
#' @rdname assoz
#' @seealso \code{\link{art1}}, \code{\link{art2}}
#' @examples 
#' \dontrun{demo(assoz_letters)}
#' \dontrun{demo(assoz_lettersSnnsR)}
#' 
#' 
#' data(snnsData)
#' patterns <- snnsData$art1_letters.pat
#' 
#' model <- assoz(patterns, dimX=7, dimY=5)
#' 
#' actMaps <- matrixToActMapList(model$fitted.values, nrow=7)
#' 
#' par(mfrow=c(3,3))
#' for (i in 1:9) plotActMap(actMaps[[i]])
assoz.default <- function(x, dimX, dimY, maxit=100, 
    initFunc="RM_Random_Weights", initFuncParams=c(1.0, -1.0), 
    learnFunc="RM_delta", learnFuncParams=c(0.01, 100, 0.0, 0.0, 0.0), 
    updateFunc="Auto_Synchronous", updateFuncParams=c(50.0),    
    shufflePatterns=TRUE, ...) {
  
  
  x <- as.matrix(x)
  
  nInputs <- dim(x)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("assoz"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)
  
  snns$archParams <- list(dimX=dimX, dimY=dimY)
  
  snns$snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Identity','Out_Identity')
  snns$snnsObject$assoz_createNet(dimX, dimY)
  
  snns <- train(snns, inputsTrain=x)
  
  #snns$fitted.values <- matrixToActMapList(snns$fitted.values, nrow=dimX)
  
  snns
}

