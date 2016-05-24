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

#' Adaptive resonance theory (ART) networks perform clustering by finding prototypes. 
#' They are mainly designed to solve the stability/plasticity dilemma (which is one of the 
#' central problems in neural networks) in the following way: new input patterns 
#' may generate new prototypes (plasticity), but patterns already present in the net 
#' (represented by their prototypes) are only altered by similar new patterns, 
#' not by others (stability).
#' ART1 is for binary inputs only,
#' if you have real-valued input, use \code{\link{art2}} instead. 
#'  
#' Learning in an ART network works as follows: 
#' A new input is intended to be classified according 
#' to the prototypes already present in the net. The similarity between the input and 
#' all prototypes is calculated. The most similar prototype is the \emph{winner}. 
#' If the similarity between the input and the winner is high enough (defined by a
#' \emph{vigilance parameter}), the winner is adapted to make it more similar to the input. 
#' If similarity is not high enough, a new prototype is created. So, at most the winner 
#' is adapted, all other prototypes remain unchanged.
#' 
#' The architecture of an ART network is the following:
#' ART is based on the more general concept of \emph{competitive learning}. The networks have 
#' two fully connected layers (in both directions), the input/comparison layer and the recognition layer. 
#' They propagate activation back and forth (resonance). The units in the recognition layer have lateral
#' inhibition, so that they show a winner-takes-all behaviour, i.e., the unit that has the highest activation
#' inhibits activation of other units, so that after a few cycles its activation will converge to one, whereas
#' the other units activations converge to zero. ART stabilizes this general learning mechanism by the presence
#' of some special units. For details refer to the referenced literature. 
#' 
#' The default initialization function, \code{ART1_Weights}, is the only one suitable for ART1 networks. It has 
#' two parameters, which are explained in the SNNS User Manual pp.189. A default of 1.0 for both is usually fine.
#' The only learning function suitable for ART1 is \code{ART1}. Update functions are \code{ART1_Stable} and 
#' \code{ART1_Synchronous}. The difference between the two is that the first one updates until the network is in a 
#' stable state, and the latter one only performs one update step. Both the learning function and the update functions 
#' have one parameter, the vigilance parameter.
#' 
#' In its current implementation, the network has two-dimensional input. The matrix \code{x} contains all 
#' (one dimensional) input patterns. Internally, every one of these patterns
#' is converted to a two-dimensional pattern using parameters \code{dimX} and \code{dimY}.
#' The parameter \code{f2Units} controls the number of units in the recognition layer, and therewith the maximal amount of clusters 
#' that are assumed to be present in the input patterns. 
#' 
#' A detailed description of the theory and the parameters is available from the SNNS documentation and the other referenced literature. 
#'
#' @title Create and train an art1 network
#' @references
#' Carpenter, G. A. & Grossberg, S. (1987), 'A massively parallel architecture for a self-organizing neural pattern recognition machine', Comput. Vision Graph. Image Process. 37, 54--115.
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
art1 <- function(x, ...) UseMethod("art1")


#' @param x a matrix with training inputs for the network
#' @param dimX x dimension of inputs and outputs
#' @param dimY y dimension of inputs and outputs
#' @param f2Units controls the number of clusters assumed to be present
#' @param maxit maximum of iterations to learn
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
# @S3method art1 default
#' @method art1 default
#' @seealso \code{\link{art2}}, \code{\link{artmap}}
#' @rdname art1
#' @examples 
#' \dontrun{demo(art1_letters)}
#' \dontrun{demo(art1_lettersSnnsR)}
#' 
#' 
#' data(snnsData)
#' patterns <- snnsData$art1_letters.pat
#' 
#' inputMaps <- matrixToActMapList(patterns, nrow=7)
#' par(mfrow=c(3,3))
#' for (i in 1:9) plotActMap(inputMaps[[i]])
#' 
#' model <- art1(patterns, dimX=7, dimY=5)
#' encodeClassLabels(model$fitted.values)
art1.default <- function(x, dimX, dimY, f2Units=nrow(x), maxit=100, 
    initFunc="ART1_Weights", initFuncParams=c(1.0, 1.0), 
    learnFunc="ART1", learnFuncParams=c(0.9, 0.0, 0.0), 
    updateFunc="ART1_Stable", updateFuncParams=c(0.0),    
    shufflePatterns=TRUE, ...) {
  
  x <- as.matrix(x)

  nInputs <- dim(x)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("art1"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)

  snns$archParams <- list(f2Units=f2Units, dimX=dimX, dimY=dimY)
  
  #snns$snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')
  snns$snnsObject$art1_createNet(dimX*dimY,dimX,f2Units,dimX)
 
  snns <- train(snns, inputsTrain=x)
  
  #snns$fitted.values <- matrixToActMapList(snns$fitted.values, nrow=dimX)
  
  snns
}

