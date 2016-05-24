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

#' Dynamic learning vector quantization (DLVQ) networks are similar to 
#' self-organizing maps (SOM, \code{\link{som}}). But they perform supervised learning
#' and lack a neighborhood relationship between the prototypes. 
#' 
#' The input data has to be normalized in order to use DLVQ.
#' 
#' Learning in DLVQ: For each class, a mean vector (prototype) is calculated and stored 
#' in a (newly generated) hidden unit. Then, the net is used to classify every pattern 
#' by using the nearest prototype. If a pattern gets misclassified as class y instead of 
#' class x, the prototype of class y is moved away from the pattern, and the prototype 
#' of class x is moved towards the pattern. This procedure is repeated iteratively until no more changes 
#' in classification take place. Then, new prototypes are introduced in the net per class
#' as new hidden units, and initialized by the mean vector of misclassified patterns in that class.
#'  
#' Network architecture: The network only has one hidden layer, containing one unit for each prototype.
#' The prototypes/hidden units are also called codebook vectors. Because SNNS generates the units 
#' automatically, and does not need their number to be specified in advance, the procedure is called
#' \emph{dynamic} LVQ in SNNS.
#' 
#' The default initialization, learning, and update functions are the only ones suitable for this kind of 
#' network. The three parameters of the learning function specify two learning rates (for the cases 
#' correctly/uncorrectly classified), and the number of cycles the net is trained before mean vectors are 
#' calculated.
#' 
#' A detailed description of the theory and the parameters is available, as always, from the SNNS 
#' documentation and the other referenced literature.
#' 
#' @title Create and train a dlvq network
#' @references 
#' Kohonen, T. (1988), Self-organization and associative memory, Vol. 8, Springer-Verlag.
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
dlvq <- function(x, ...) UseMethod("dlvq")


#' @param x a matrix with training inputs for the network
#' @param y the corresponding target values
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
# @S3method dlvq default
#' @method dlvq default
#' @rdname dlvq
#' @examples 
#' \dontrun{demo(dlvq_ziff)}
#' \dontrun{demo(dlvq_ziffSnnsR)}
#' 
#' 
#' data(snnsData)
#' dataset <- snnsData$dlvq_ziff_100.pat
#' 
#' inputs <- dataset[,inputColumns(dataset)]
#' outputs <- dataset[,outputColumns(dataset)]
#' 
#' model <- dlvq(inputs, outputs)
#' 
#' fitted(model) == outputs
#' mean(fitted(model) - outputs)
dlvq.default <- function(x, y, 
    initFunc="DLVQ_Weights", initFuncParams=c(1.0, -1.0), 
    learnFunc="Dynamic_LVQ", learnFuncParams=c(0.03, 0.03, 10.0), 
    updateFunc="Dynamic_LVQ", updateFuncParams=c(0.0),    
    shufflePatterns=TRUE, ...) {
  
  
  x <- as.matrix(x)
  
  nInputs <- dim(x)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("dlvq"), nInputs=nInputs, maxit=1, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)
  
  snns$archParams <- NULL #list(dimX=dimX, dimY=dimY)
  
  snns$snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Identity','Out_Identity')
  snns$snnsObject$createNet(c(nInputs, 1), fullyConnectedFeedForward = FALSE)
  
  snns <- train(snns, inputsTrain=x, targetsTrain=y)
  
  #snns$fitted.values <- matrixToActMapList(snns$fitted.values, nrow=dimX)
  
  snns
}

