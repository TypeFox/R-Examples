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


#' Print out some characteristics of an \code{\link{rsnns}} object.
#' 
#' @title Generic print function for rsnns objects
#' @param x the \code{\link{rsnns}} object
#' @param ... additional function parameters (currently not used)
#' @export
# @S3method print rsnns
#' @method print rsnns
# @rdname rsnns
print.rsnns <- function(x, ...) {
  if(!inherits(x, "rsnns")) stop("not a legitimate rsnns model")
  
  cat("Class: ", paste(class(x), sep="", collapse="->"), "\n", sep="")
  cat("Number of inputs:",x$nInputs, "\n",sep=" ")
  cat("Number of outputs:",x$nOutputs, "\n",sep=" ")
  cat("Maximal iterations:",x$maxit, "\n",sep=" ")
  cat("Initialization function:",x$initFunc, "\n",sep=" ")
  cat("Initialization function parameters:",x$initFuncParams, "\n",sep=" ")
  cat("Learning function:",x$learnFunc, "\n",sep=" ")  
  cat("Learning function parameters:",x$learnFuncParams, "\n",sep=" ")
  cat("Update function:",x$updateFunc, "\n",sep="")
  cat("Update function parameters:",x$updateFuncParams, "\n",sep=" ")  
  cat("Patterns are shuffled internally:",x$shufflePatterns, "\n",sep=" ")
  cat("Compute error in every iteration:",x$computeIterativeError, "\n",sep=" ")  
  cat("Architecture Parameters:\n",sep="")
  print(x$archParams)
  cat("All members of model:\n",sep="") 
  print(names(x))
  
  invisible(x)
}


#' This function generates a list of data.frames containing the most important information 
#' that defines a network, in a format that is easy to use. To get the full definition in 
#' the original SNNS format, use \code{\link{summary.rsnns}} or \code{\link{exportToSnnsNetFile}} 
#' instead. 
#' 
#' Internally, a call to \code{\link{SnnsRObject$extractNetInfo}} is done, and the results of 
#' this call are returned.
#' 
#' @title Extract information from a network
#' @param object the \code{\link{rsnns}} object
#' @return a list containing information extracted from the network (see \code{\link{SnnsRObject$extractNetInfo}}).
#' @export
#' @seealso \code{\link{SnnsRObject$extractNetInfo}}
extractNetInfo <- function(object) {
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  object$snnsObject$extractNetInfo()
}

#' Export the net that is present in the \code{\link{rsnns}} object in the 
#' original (.net) SNNS file format.
#'
#' @title Export the net to a file in the original SNNS file format
#' @param object the \code{\link{rsnns}} object
#' @param filename path and filename to be written to 
#' @param netname name that is given to the network in the file
#' @export
exportToSnnsNetFile <- function(object, filename, netname="RSNNS_untitled") {
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  object$snnsObject$saveNet(filename, netname)
}


#' Prints out a summary of the network. The printed information can be either 
#' all information of the network in the original SNNS file format,
#' or the information given by \code{\link{extractNetInfo}}.
#' This behaviour is controlled with the parameter \code{origSnnsFormat}.
#' 
#' @title Generic summary function for rsnns objects
#' @param object the \code{\link{rsnns}} object
#' @param origSnnsFormat show data in SNNS's original format in which networks are saved, or show output of \code{\link{extractNetInfo}}
#' @param ... additional function parameters (currently not used)
#' @return Either the contents of the .net file that SNNS would generate from 
#' the object, as a string. Or the output of \code{\link{extractNetInfo}}.  
#' @export
# @S3method summary rsnns
#' @method summary rsnns
#' @seealso \code{\link{extractNetInfo}} 
summary.rsnns <- function(object, origSnnsFormat=TRUE, ...) {
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  if(origSnnsFormat) {

    s <- object$snnsObject$serializeNet("RSNNS_untitled")
    s <- s$serialization
#    filename <- tempfile(pattern = "rsnns")
#    object$snnsObject$saveNet(filename, " ")
#    file <- file(filename, "r")
#    s <- readLines(file)
#    close(file)
#    unlink(filename)    
  } else {
    s <- extractNetInfo(object)
    
    if(length(s$fullWeightMatrix) > 20*20)
      s$fullWeightMatrix <- "omitting full weight matrix as it is bigger than 20*20"
      
  }
  #invisible(s$serialization)
  cat(s)
  invisible(s)
}


#' The object factory generates an \code{rsnns} object and initializes its
#' member variables with the values given as parameters. Furthermore, it
#' generates an object of \code{\link{SnnsR-class}}. Later, this information is
#' to be used to train the network.
#' 
#' The typical procedure implemented in \code{rsnns} subclasses is the following: 
#' \itemize{
#' \item generate the \code{rsnns} object with this object factory
#' \item generate the network according to the architecture needed
#' \item train the network (with \code{\link{train}})
#' }
#'
#' In every \code{rsnns} object, the iterative error is the summed squared error
#' (SSE) of all patterns. If the SSE is computed on the test set, then it is
#' weighted to take care of the different amount of patterns in the sets.
#' 
#' @title Object factory for generating rsnns objects
#' @param subclass the subclass of rsnns to generate (vector of strings)
#' @param nInputs the number of inputs the network will have
#' @param maxit maximum of iterations to learn
#' @param initFunc the initialization function to use
#' @param initFuncParams the parameters for the initialization function
#' @param learnFunc the learning function to use
#' @param learnFuncParams the parameters for the learning function
#' @param updateFunc the update function to use
#' @param updateFuncParams the parameters for the update function
#' @param shufflePatterns should the patterns be shuffled?
#' @param computeIterativeError should the error be computed in every iteration? 
#' @param pruneFunc the pruning function to use
#' @param pruneFuncParams the parameters for the pruning function. Unlike the other functions, 
#' these have to be given in a named list. See the pruning demos for further explanation. 
#' @return a partly initialized \code{rsnns} object 
#' @aliases rsnns
#' @export
#' @seealso \code{\link{mlp}}, \code{\link{dlvq}}, \code{\link{rbf}}, \code{\link{rbfDDA}}, \code{\link{elman}}, 
#' \code{\link{jordan}}, \code{\link{som}}, \code{\link{art1}}, \code{\link{art2}}, \code{\link{artmap}}, \code{\link{assoz}}
rsnnsObjectFactory <- function(subclass, nInputs, maxit, 
    initFunc, initFuncParams, 
    learnFunc, learnFuncParams, 
    updateFunc, updateFuncParams, 
    shufflePatterns=TRUE, computeIterativeError=TRUE, pruneFunc=NULL, pruneFuncParams=NULL) {
  
  snns <- NULL
  
  snns$nInputs <- nInputs
  snns$maxit <- maxit  
  
  snns$initFunc <- initFunc
  snns$initFuncParams <- initFuncParams 
  snns$learnFunc <- learnFunc
  snns$learnFuncParams <- learnFuncParams  
  snns$updateFunc <- updateFunc
  snns$updateFuncParams <- updateFuncParams    
  snns$shufflePatterns <- shufflePatterns
  snns$computeIterativeError <- computeIterativeError
  
  snns$pruneFunc <- pruneFunc
  snns$pruneFuncParams <- pruneFuncParams
  
  snns$snnsObject <- SnnsRObjectFactory()
  
  class(snns) <- c(subclass, "rsnns")
    
  snns
}


#' The function calls \code{\link{SnnsRObject$train}} and saves the result in the
#' current \code{\link{rsnns}} object. This function is used internally by the 
#' models (e.g. \code{\link{mlp}}) for training. Unless you are not about to implement
#' a new model on the S3 layer you most probably don't want to use this function.
#' 
#' @title Internal generic train function for rsnns objects
# @param object the object to which to apply train
# @param ... additional function parameters
#' @export
train <- function(object, ...) UseMethod("train")

#' Internal generic train function for \code{rsnns} objects.
#' 
#' @param object the \code{\link{rsnns}} object
#' @param inputsTrain training input
#' @param targetsTrain training targets
#' @param inputsTest test input
#' @param targetsTest test targets
#' @param serializeTrainedObject parameter passed to \code{\link{SnnsRObject$train}}
#' @param ... additional function parameters (currently not used)
#' @return an \code{\link{rsnns}} object, to which the results of training have been added. 
#' @export
# @S3method train rsnns
#' @method train rsnns
#' @rdname train
train.rsnns <- function(object, inputsTrain, targetsTrain=NULL, inputsTest=NULL, targetsTest=NULL, serializeTrainedObject=TRUE, ...) {
  
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  trainResult <- object$snnsObject$train(inputsTrain=inputsTrain, targetsTrain=targetsTrain, 
      initFunc=object$initFunc, initFuncParams=object$initFuncParams, 
      learnFunc=object$learnFunc, learnFuncParams=object$learnFuncParams, updateFunc=object$updateFunc, 
      updateFuncParams=object$updateFuncParams, outputMethod=class(object)[1], maxit=object$maxit, 
      shufflePatterns=object$shufflePatterns, computeError=object$computeIterativeError, 
      inputsTest=inputsTest, targetsTest=targetsTest, serializeTrainedObject=serializeTrainedObject,
      pruneFunc=object$pruneFunc, pruneFuncParams=object$pruneFuncParams)
  
  object$IterativeFitError <- trainResult$IterativeFitError
  object$IterativeTestError <- trainResult$IterativeTestError
  
  object$fitted.values <- trainResult$fitValues
  object$fittedTestValues <- trainResult$testValues
  
  object$nOutputs <- ncol(trainResult$fitValues)
  object
}


#' Predict values using the given network. 
#'
#' @title Generic predict function for rsnns object
#' @param object the \code{\link{rsnns}} object
#' @param newdata the new input data which is used for prediction
#' @param ... additional function parameters (currently not used)
#' @return the predicted values
# @S3method predict rsnns
#' @method predict rsnns
# @rdname rsnns
#' @export
predict.rsnns <- function(object, newdata, ...) {
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  #type <- match.arg(type)
  if(missing(newdata)) z <- fitted(object)
  else {
    if(is.null(dim(newdata)))
      dim(newdata) <- c(1L, length(newdata)) # a row vector
    x <- as.matrix(newdata)     # to cope with dataframes
    if(any(is.na(x))) stop("missing values in 'x'")
    
    keep <- 1L:nrow(x)
    rn <- rownames(x)
    
    ntr <- nrow(x)
    nout <- object$nOutputs
    
    z <- matrix(NA, nrow(newdata), nout, dimnames = list(rn, dimnames(object$fitted.values)[[2L]]))
    
    #predict values.. 
    patSet <- object$snnsObject$createPatSet(newdata) 
    predictions <- object$snnsObject$predictCurrPatSet(class(object)[1], updateFuncParams=object$updateFuncParams)
    object$snnsObject$deletePatSet(patSet$set_no)
    z[keep,] <- predictions
  }
  z
  #predictions
}


#' The function calls \code{\link{SnnsRObject$getCompleteWeightMatrix}} and returns its result.
#'
#' @title Function to extract the weight matrix of an rsnns object
# @param object the object to which to apply weightMatrix
# @param ... additional function parameters
#' @export
weightMatrix <- function(object, ...) UseMethod("weightMatrix")

#' Function to extract the weight matrix of an rsnns object.
#' 
#' @param object the \code{\link{rsnns}} object
#' @param ... additional function parameters (currently not used)
#' @return a matrix with all weights from all neurons present in the net. 
#' @export
# @S3method weightMatrix rsnns
#' @method weightMatrix rsnns
#' @rdname weightMatrix
weightMatrix.rsnns <- function(object, ...) {
  
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  object$snnsObject$getCompleteWeightMatrix(setDimNames=TRUE)
}
