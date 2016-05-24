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


#' Jordan networks are partially recurrent networks and similar to Elman
#' networks (see \code{\link{elman}}). Partially recurrent networks are useful
#' when working with time series data. I.e., when the output of the network not
#' only should depend on the current pattern, but also on the patterns presented
#' before.
#' 
#' Learning on Jordan networks: Backpropagation algorithms for feed-forward
#' networks can be adapted for their use with this type of networks. In SNNS,
#' there exist adapted versions of several backpropagation-type algorithms for
#' Jordan and Elman networks.
#' 
#' Network architecture: A Jordan network can be seen as a feed-forward network
#' with additional context units in the input layer. These context units take
#' input from themselves (direct feedback), and from the output units. The
#' context units save the current state of the net. In a Jordan net, the number
#' of context units and output units has to be the same.
#' 
#' Initialization of Jordan and Elman nets should be done with the default init
#' function \code{JE_Weights}, which has five parameters. The first two
#' parameters define an interval from which the forward connections are randomly
#' chosen. The third parameter gives the self-excitation weights of the context
#' units. The fourth parameter gives the weights of context units between them,
#' and the fifth parameter gives the initial activation of context units.
#' 
#' Learning functions are \code{JE_BP}, \code{JE_BP_Momentum},
#' \code{JE_Quickprop}, and \code{JE_Rprop}, which are all adapted versions of
#' their standard-procedure counterparts.  Update functions that can be used are
#' \code{JE_Order} and \code{JE_Special}.
#' 
#' A detailed description of the theory and the parameters is available, as
#' always, from the SNNS documentation and the other referenced literature.
#' 
#' @title Create and train a Jordan network
#' @references Jordan, M. I. (1986), 'Serial Order: A Parallel, Distributed
#' Processing Approach', Advances in Connectionist Theory Speech 121(ICS-8604),
#' 471-495.
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual,
#' Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen.
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
jordan <- function(x, ...) UseMethod("jordan")


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
#' @param inputsTest a matrix with inputs to test the network
#' @param targetsTest the corresponding targets for the test input
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object.
#' @export
# @S3method jordan default
#' @method jordan default
#' @rdname jordan
#' @seealso \code{\link{elman}}
#' @examples 
#' \dontrun{demo(iris)}
#' \dontrun{demo(laser)}
#' \dontrun{demo(eight_elman)}
#' \dontrun{demo(eight_elmanSnnsR)}
#' 
#' 
#' data(snnsData)
#' inputs <- snnsData$laser_1000.pat[,inputColumns(snnsData$laser_1000.pat)]
#' outputs <- snnsData$laser_1000.pat[,outputColumns(snnsData$laser_1000.pat)]
#' 
#' patterns <- splitForTrainingAndTest(inputs, outputs, ratio=0.15)
#' 
#' modelJordan <- jordan(patterns$inputsTrain, patterns$targetsTrain, 
#'                        size=c(8), learnFuncParams=c(0.1), maxit=100,
#'                        inputsTest=patterns$inputsTest, 
#'                        targetsTest=patterns$targetsTest, linOut=FALSE)
#' 
#' names(modelJordan)
#' 
#' par(mfrow=c(3,3))
#' plotIterativeError(modelJordan)
#' 
#' plotRegressionError(patterns$targetsTrain, modelJordan$fitted.values)
#' plotRegressionError(patterns$targetsTest, modelJordan$fittedTestValues)
#' hist(modelJordan$fitted.values - patterns$targetsTrain, col="lightblue")
#' 
#' plot(inputs, type="l")
#' plot(inputs[1:100], type="l")
#' lines(outputs[1:100], col="red")
#' lines(modelJordan$fitted.values[1:100], col="green")
jordan.default <- function(x, y, size=c(5), maxit=100, 
    initFunc="JE_Weights", initFuncParams=c(1.0,  -1.0,  0.3,  1.0,  0.5), 
    learnFunc="JE_BP", learnFuncParams=c(0.2), 
    updateFunc="JE_Order", updateFuncParams=c(0.0),    
    shufflePatterns=FALSE, linOut=TRUE, inputsTest=NULL, targetsTest=NULL, ...) {
  
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  checkInput(x,y)
  
  nInputs <- dim(x)[2L]
  nOutputs <- dim(y)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("jordan"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=TRUE)
  
  snns$archParams <- list(size=size)
  
  snns$snnsObject$setUnitDefaults(1,0,1,0,1,"Act_Logistic","Out_Identity")
  snns$snnsObject$jordan_createNet(nInputs, size, nOutputs, 1, 1, 1)
  
  
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
