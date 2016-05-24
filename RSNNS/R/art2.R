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


#' ART2 is very similar to ART1, but for real-valued input. See \code{\link{art1}}
#' for more information. Opposed to the ART1 implementation, the ART2 implementation 
#' does not assume two-dimensional input. 
#'
#' As comparison of real-valued vectors is more difficult than comparison of binary vectors, 
#' the comparison layer is more complex in ART2, and actually consists of three layers. With a more
#' complex comparison layer, also other parts of the network enhance their complexity.
#' In SNNS, this enhanced complexity is reflected by the presence of more parameters in initialization-, learning-,
#' and update function.  
#' 
#' In analogy to the implementation of ART1, there are one initialization function, one learning function and two 
#' update functions suitable for ART2.  The learning and update functions have five parameters, the initialization function has two 
#' parameters. For details see the SNNS User Manual, p. 67 and pp. 192.
#'  
#' 
#' @title Create and train an art2 network
#' @references
#' Carpenter, G. A. & Grossberg, S. (1987), 'ART 2: self-organization of stable category recognition codes for analog input patterns', Appl. Opt. 26(23), 4919--4930.
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
art2 <- function(x, ...) UseMethod("art2")


#' @param x a matrix with training inputs for the network
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
#' @return an \code{\link{rsnns}} object. The \code{fitted.values} member contains the 
#' activation patterns for all inputs.
#' @export
# @S3method art2 default
#' @method art2 default
#' @seealso \code{\link{art1}}, \code{\link{artmap}}
#' @rdname art2
#' @examples 
#' \dontrun{demo(art2_tetra)}
#' \dontrun{demo(art2_tetraSnnsR)}
#' 
#' 
#' data(snnsData)
#' patterns <- snnsData$art2_tetra_med.pat
#' 
#' model <- art2(patterns, f2Units=5, learnFuncParams=c(0.99, 20, 20, 0.1, 0), 
#'                   updateFuncParams=c(0.99, 20, 20, 0.1, 0))
#' model
#' 
#' testPatterns <- snnsData$art2_tetra_high.pat
#' predictions <- predict(model, testPatterns)
#' 
#' \dontrun{library(scatterplot3d)}
#' 
#' \dontrun{par(mfrow=c(2,2))}
#' \dontrun{scatterplot3d(patterns, pch=encodeClassLabels(model$fitted.values))}
#' \dontrun{scatterplot3d(testPatterns, pch=encodeClassLabels(predictions))}
art2.default <- function(x, f2Units=5, maxit=100, 
    initFunc="ART2_Weights", initFuncParams=c(0.9, 2.0), 
    learnFunc="ART2", learnFuncParams=c(0.98, 10.0, 10.0, 0.1, 0.0), 
    updateFunc="ART2_Stable", updateFuncParams=c(0.98, 10.0, 10.0, 0.1, 0.0),    
    shufflePatterns=TRUE, ...) {
  
  
  x <- as.matrix(x)
  
  nInputs <- dim(x)[2L]
  
  snns <- rsnnsObjectFactory(subclass=c("art2"), nInputs=nInputs, maxit=maxit, 
      initFunc=initFunc, initFuncParams=initFuncParams, 
      learnFunc=learnFunc, learnFuncParams=learnFuncParams, 
      updateFunc=updateFunc, 
      updateFuncParams=updateFuncParams,
      shufflePatterns=shufflePatterns, computeIterativeError=FALSE)
  
  snns$archParams <- list(f2Units=f2Units)
  
  #snns$snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')
  snns$snnsObject$art2_createNet(nInputs,nInputs,f2Units,f2Units)

  snns <- train(snns, inputsTrain=x)
  
  snns
}

