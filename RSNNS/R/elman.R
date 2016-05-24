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


#' Elman networks are partially recurrent networks and similar to Jordan
#' networks (function \code{\link{jordan}}). For details, see explanations
#' there.
#' 
#' Learning in Elman networks:
#' Same as in Jordan networks (see \code{\link{jordan}}).
#' 
#' Network architecture: The difference between Elman and Jordan networks is
#' that in an Elman network the context units get input not from the output
#' units, but from the hidden units. Furthermore, there is no direct feedback in
#' the context units. In an Elman net, the number of context units and hidden
#' units has to be the same. The main advantage of Elman nets is that the number
#' of context units is not directly determined by the output dimension (as in
#' Jordan nets), but by the number of hidden units, which is more flexible, as
#' it is easy to add/remove hidden units, but not output units.
#' 
#' A detailed description of the theory and the parameters is available, as
#' always, from the SNNS documentation and the other referenced literature.
#' 
#' @title Create and train an Elman network
#' @references
#'
#' Elman, J. L. (1990), 'Finding structure in time', Cognitive Science 14(2),
#' 179--211.
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual,
#' Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen.
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
elman <- function(x, ...) UseMethod("elman")


#' @param x a matrix with training inputs for the network
#' @param y the corresponding targets values
#' @param size number of units in the hidden layer(s)
#' @param maxit maximum of iterations to learn
#' @param initFunc the initialization function to use
#' @param initFuncParams the parameters for the initialization function
#' @param learnFunc the learning function to use
#' @param learnFuncParams the parameters for the learning function
#' @param updateFunc the update function to use
#' @param updateFuncParams the parameters for the update function
#' @param shufflePatterns should the patterns be shuffled?
#' @param linOut sets the activation function of the output units to linear or logistic
#' @param outContext if TRUE, the context units are also output units (untested)
#' @param inputsTest a matrix with inputs to test the network
#' @param targetsTest the corresponding targets for the test input
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object.
#' @export
# @S3method elman default
#' @method elman default
#' @rdname elman
#' @seealso \code{\link{jordan}}
#' @examples 
#' \dontrun{demo(iris)}
#' \dontrun{demo(laser)}
#' \dontrun{demo(eight_elman)}
#' \dontrun{demo(eight_elmanSnnsR)}
#' 
#' 
#' data(snnsData)
#' inputs <- snnsData$eight_016.pat[,inputColumns(snnsData$eight_016.pat)]
#' outputs <- snnsData$eight_016.pat[,outputColumns(snnsData$eight_016.pat)]
#' 
#' par(mfrow=c(1,2))
#' 
#' modelElman <- elman(inputs, outputs, size=8, learnFuncParams=c(0.1), maxit=1000)
#' modelElman
#' modelJordan <- jordan(inputs, outputs, size=8, learnFuncParams=c(0.1), maxit=1000)
#' modelJordan
#' 
#' plotIterativeError(modelElman)
#' plotIterativeError(modelJordan)
#' 
#' summary(modelElman)
#' summary(modelJordan)
elman.default <- function(x, y, size=c(5), maxit=100, 
    initFunc="JE_Weights", initFuncParams=c(1.0,  -1.0,  0.3,  1.0,  0.5), 
    learnFunc="JE_BP", learnFuncParams=c(0.2), 
    updateFunc="JE_Order", updateFuncParams=c(0.0),    
    shufflePatterns=FALSE, linOut=TRUE, outContext=FALSE, inputsTest=NULL, targetsTest=NULL, ...) {
  

  x <- as.matrix(x)
  y <- as.matrix(y)
  
  checkInput(x,y)
  
  nInputs <- dim(x)[2L]
  nOutputs <- dim(y)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("elman"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=TRUE)
  
  snns$archParams <- list(size=size)
  
  snns$snnsObject$setUnitDefaults(1,0,1,0,1,"Act_Logistic","Out_Identity")
  snns$snnsObject$elman_createNet(c(nInputs, size, nOutputs), seq(1,1,length=(length(size)+2)), outContext)
  
  if(linOut) {
    outputActFunc <- "Act_Identity"
  } else {
    outputActFunc <- "Act_Logistic"
  }
  
  snns$snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
  snns$snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", outputActFunc)
  
  
  snns <- train(snns, inputsTrain=x, targetsTrain=y, inputsTest=inputsTest, targetsTest=targetsTest)

  snns  
}
