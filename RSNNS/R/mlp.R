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


#' This function creates a multilayer perceptron (MLP) and trains it. MLPs are
#' fully connected feedforward networks, and probably the most common network
#' architecture in use.  Training is usually performed by error backpropagation
#' or a related procedure.
#'  
#' There are a lot of different learning functions present in SNNS that can be
#' used together with this function, e.g., \code{Std_Backpropagation},
#' \code{BackpropBatch}, \code{BackpropChunk}, \code{BackpropMomentum},
#' \code{BackpropWeightDecay}, \code{Rprop}, \code{Quickprop}, \code{SCG}
#' (scaled conjugate gradient), ...
#' 
#' \code{Std_Backpropagation}, \code{BackpropBatch}, e.g., have two parameters,
#' the learning rate and the maximum output difference. The learning rate is
#' usually a value between 0.1 and 1. It specifies the gradient descent step
#' width. The maximum difference defines, how much difference between output and
#' target value is treated as zero error, and not backpropagated. This parameter
#' is used to prevent overtraining. For a complete list of the parameters of all
#' the learning functions, see the SNNS User Manual, pp. 67.
#' 
#' The defaults that are set for initialization and update functions usually don't have to be changed. 
#' 
#' 
#' 
#' @title Create and train a multi-layer perceptron (MLP)
#' @references Rosenblatt, F. (1958), 'The perceptron: A probabilistic model for
#' information storage and organization in the brain', Psychological Review
#' 65(6), 386--408.
#' 
#' Rumelhart, D. E.; Clelland, J. L. M. & Group, P. R. (1986), Parallel distributed processing :explorations in the microstructure of cognition, Mit, Cambridge, MA etc.
#'  
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
mlp <- function(x, ...) UseMethod("mlp")


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
#' @param hiddenActFunc the activation function of all hidden units
#' @param shufflePatterns should the patterns be shuffled?
#' @param linOut sets the activation function of the output units to linear or logistic
#' @param inputsTest a matrix with inputs to test the network
#' @param targetsTest the corresponding targets for the test input
#' @param pruneFunc the pruning function to use
#' @param pruneFuncParams the parameters for the pruning function. Unlike the other functions, 
#' these have to be given in a named list. See the pruning demos for further explanation. 
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object.
#' @export
# @S3method mlp default
#' @method mlp default
#' @rdname mlp
#' @examples 
#' \dontrun{demo(iris)}
#' \dontrun{demo(laser)}
#' \dontrun{demo(encoderSnnsCLib)}
#' 
#' 
#' data(iris)
#' 
#' #shuffle the vector
#' iris <- iris[sample(1:nrow(iris),length(1:nrow(iris))),1:ncol(iris)]
#' 
#' irisValues <- iris[,1:4]
#' irisTargets <- decodeClassLabels(iris[,5])
#' #irisTargets <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)
#' 
#' iris <- splitForTrainingAndTest(irisValues, irisTargets, ratio=0.15)
#' iris <- normTrainingAndTestSet(iris)
#' 
#' model <- mlp(iris$inputsTrain, iris$targetsTrain, size=5, learnFuncParams=c(0.1), 
#'               maxit=50, inputsTest=iris$inputsTest, targetsTest=iris$targetsTest)
#' 
#' summary(model)
#' model
#' weightMatrix(model)
#' extractNetInfo(model)
#' 
#' par(mfrow=c(2,2))
#' plotIterativeError(model)
#' 
#' predictions <- predict(model,iris$inputsTest)
#' 
#' plotRegressionError(predictions[,2], iris$targetsTest[,2])
#' 
#' confusionMatrix(iris$targetsTrain,fitted.values(model))
#' confusionMatrix(iris$targetsTest,predictions)
#' 
#' plotROC(fitted.values(model)[,2], iris$targetsTrain[,2])
#' plotROC(predictions[,2], iris$targetsTest[,2])
#' 
#' #confusion matrix with 402040-method
#' confusionMatrix(iris$targetsTrain, encodeClassLabels(fitted.values(model),
#'                                                        method="402040", l=0.4, h=0.6))
mlp.default <- function(x, y, size=c(5), maxit=100,  
    initFunc="Randomize_Weights", initFuncParams=c(-0.3, 0.3), 
    learnFunc="Std_Backpropagation", learnFuncParams=c(0.2, 0.0), 
    updateFunc="Topological_Order", updateFuncParams=c(0.0),
    hiddenActFunc="Act_Logistic",
    shufflePatterns=TRUE, linOut=FALSE, inputsTest=NULL, targetsTest=NULL, pruneFunc=NULL, pruneFuncParams=NULL, ...) {


  x <- as.matrix(x)
  y <- as.matrix(y)

  checkInput(x, y)

  nInputs <- ncol(x)
  nOutputs <- ncol(y)


  snns <- rsnnsObjectFactory(
    subclass=c("mlp"), nInputs=nInputs, maxit=maxit, 
    initFunc=initFunc, initFuncParams=initFuncParams, 
    learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
    updateFunc=updateFunc, updateFuncParams=updateFuncParams,
    shufflePatterns=shufflePatterns, computeIterativeError=TRUE,
    pruneFunc=pruneFunc, pruneFuncParams=pruneFuncParams)
  
  snns$archParams <- list(size=size)
  snns$snnsObject$setUnitDefaults(0,0,1,0,1,"Act_Logistic","Out_Identity")
  snns$snnsObject$createNet(unitsPerLayer=c(nInputs, size, nOutputs),
                            fullyConnectedFeedForward=TRUE,
                            iNames = colnames(x), oNames = colnames(y))
  
  outputActFunc <- if(linOut)  "Act_Identity" else  "Act_Logistic"
  
  snns$snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
  snns$snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", hiddenActFunc)
  snns$snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", outputActFunc)
  
  snns <- train(snns,
                inputsTrain = x, inputsTest = inputsTest,
                targetsTrain = y, targetsTest = targetsTest)
  
}

