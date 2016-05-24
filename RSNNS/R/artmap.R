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


#' An ARTMAP performs supervised learning. It consists of two coupled ART networks.
#' In theory, these could be ART1, ART2, or others. However, in SNNS ARTMAP is
#' implemented for ART1 only. So, this function is to be used with binary input. 
#' As explained in the description of \code{\link{art1}}, ART aims at solving the stability/plasticity
#' dilemma. So the advantage of ARTMAP is that it is a supervised learning mechanism
#' that guarantees stability.
#'
#' See also the details section of \code{\link{art1}}. The two ART1 networks are connected by a \emph{map field}.
#' The input of the first ART1 network is the training input, the input of the second network are the target values, 
#' the teacher signals. The two networks are often called ARTa and ARTb, we call them here training data network 
#' and target data network.  
#' 
#' In analogy to the ART1 and ART2 implementations, there are one initialization function, one learning function, 
#' and two update functions present that are suitable for ARTMAP. The parameters are basically as in ART1, but for 
#' two networks. The learning function and the update functions have 3 parameters, the vigilance parameters of the 
#' two ART1 networks and an additional vigilance parameter for inter ART reset control. The initialization function
#' has four parameters, two for every ART1 network.
#' 
#' A detailed description of the theory and the parameters is available from 
#' the SNNS documentation and the other referenced literature. 
#' 
#' @title Create and train an artmap network
#' @references
#' Carpenter, G. A.; Grossberg, S. & Reynolds, J. H. (1991), 'ARTMAP: Supervised real-time learning and classification of nonstationary data by a self-organizing neural network', Neural Networks 4(5), 565--588.
#' 
#' Grossberg, S. (1988), Adaptive pattern classification and universal recoding. I.: parallel development and coding of neural feature detectors, MIT Press, Cambridge, MA, USA, chapter I, pp. 243--258.
#' 
#' Herrmann, K.-U. (1992), 'ART -- Adaptive Resonance Theory -- Architekturen, Implementierung und Anwendung', Master's thesis, IPVR, University of Stuttgart. (in German)
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley. (in German)
#' @export
artmap <- function(x, ...) UseMethod("artmap")



#' @param x a matrix with training inputs and targets for the network
#' @param nInputsTrain the number of columns of the matrix that are training input
#' @param nInputsTargets the number of columns that are target values
#' @param nUnitsRecLayerTrain number of units in the recognition layer of the training data ART network
#' @param nUnitsRecLayerTargets number of units in the recognition layer of the target data ART network
#' @param maxit maximum of iterations to perform
#' @param nRowInputsTrain number of rows the training input units are to be organized in (only for visualization purposes of the net in the original SNNS software)
#' @param nRowInputsTargets same, but for the target value input units
#' @param nRowUnitsRecLayerTrain same, but for the recognition layer of the training data ART network
#' @param nRowUnitsRecLayerTargets same, but for the recognition layer of the target data ART network
#' @param initFunc the initialization function to use
#' @param initFuncParams the parameters for the initialization function
#' @param learnFunc the learning function to use
#' @param learnFuncParams the parameters for the learning function
#' @param updateFunc the update function to use
#' @param updateFuncParams the parameters for the update function
#' @param shufflePatterns should the patterns be shuffled?
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object. The \code{fitted.values} member of the object contains a 
#' list of two-dimensional activation patterns.
#' @export
# @S3method artmap default
#' @method artmap default
#' @seealso \code{\link{art1}}, \code{\link{art2}}
#' @rdname artmap
#' @examples 
#' \dontrun{demo(artmap_letters)}
#' \dontrun{demo(artmap_lettersSnnsR)}
#' 
#' 
#' data(snnsData)
#' trainData <- snnsData$artmap_train.pat
#' testData <- snnsData$artmap_test.pat
#' 
#' model <- artmap(trainData, nInputsTrain=70, nInputsTargets=5, 
#'                   nUnitsRecLayerTrain=50, nUnitsRecLayerTargets=26)
#' model$fitted.values
#' 
#' predict(model, testData)
artmap.default <- function(x, nInputsTrain, nInputsTargets, nUnitsRecLayerTrain, nUnitsRecLayerTargets, maxit=1, 
    nRowInputsTrain=1, nRowInputsTargets=1, nRowUnitsRecLayerTrain=1, nRowUnitsRecLayerTargets=1,
    initFunc="ARTMAP_Weights", initFuncParams=c(1.0, 1.0, 1.0, 1.0, 0.0), 
    learnFunc="ARTMAP", learnFuncParams=c(0.8, 1.0, 1.0, 0, 0), 
    updateFunc="ARTMAP_Stable", updateFuncParams=c(0.8, 1.0, 1.0, 0, 0),    
    shufflePatterns=TRUE, ...) {
  
  x <- as.matrix(x)

  nInputs <- nInputsTrain + nInputsTargets
  
  if(nInputs != dim(x)[2L]) {
    warning("nInputsTrain + nInputsTargets has to be the same as the number of columns of the input matrix. Not training the ARTMAP..")
    return()
  }

  snns <- rsnnsObjectFactory(subclass=c("artmap"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)

  snns$archParams <- list(nInputsTrain=nInputsTrain, nInputsTargets=nInputsTargets, 
                          nUnitsRecLayerTrain=nUnitsRecLayerTrain, nUnitsRecLayerTargets=nUnitsRecLayerTargets, 
      nRowInputsTrain=nRowInputsTrain, nRowInputsTargets=nRowInputsTargets, 
      nRowUnitsRecLayerTrain=nRowUnitsRecLayerTrain, nRowUnitsRecLayerTargets=nRowUnitsRecLayerTargets)
    
  snns$snnsObject$artmap_createNet(nInputsTrain,nRowInputsTrain,nUnitsRecLayerTrain,nRowUnitsRecLayerTrain,
                              nInputsTargets,nRowInputsTargets,nUnitsRecLayerTargets,nRowUnitsRecLayerTargets)
  
  snns <- train(snns, inputsTrain=x)
  
  snns
}

