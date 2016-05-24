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



#' Split the input and target values to a training and a test set. Test set is taken from the end of the
#' data. If the data is to be shuffled, this should be done before calling this function.
#' 
#' @title Function to split data into training and test set
#' @param x inputs
#' @param y targets
#' @param ratio ratio of training and test sets (default: 15\% of the data is used for testing)
#' @return a named list with the following elements:
#' \item{inputsTrain}{a matrix containing the training inputs}
#' \item{targetsTrain}{a matrix containing the training targets}
#' \item{inputsTest}{a matrix containing the test inputs}
#' \item{targetsTest}{a matrix containing the test targets}
#' @export
#' @examples
#' data(iris)
#' #shuffle the vector
#' iris <- iris[sample(1:nrow(iris),length(1:nrow(iris))),1:ncol(iris)]
#' 
#' irisValues <- iris[,1:4]
#' irisTargets <- decodeClassLabels(iris[,5])
#' 
#' splitForTrainingAndTest(irisValues, irisTargets, ratio=0.15)
splitForTrainingAndTest <- function(x, y, ratio=0.15) {
  
  x <- as.matrix(x)
  nInputs <- nrow(x)
  
  y <- as.matrix(y)
  
  trainIndices <- 1 : (nInputs * (1-ratio))
  testIndices <- (1:nInputs)[-trainIndices]
  
  inputsTrain <- x[trainIndices,]
  targetsTrain <- y[trainIndices,]
  
  inputsTest <- x[testIndices,]
  targetsTest <- y[testIndices,]
  
  list(inputsTrain=inputsTrain, targetsTrain=targetsTrain, inputsTest=inputsTest, targetsTest=targetsTest)  
}

#' Normalize training and test set as obtained by \code{\link{splitForTrainingAndTest}} in the following way:
#' The \code{inputsTrain} member is normalized using \code{\link{normalizeData}} with the parameters given in \code{type}.
#' The normalization parameters obtained during this normalization are then used to normalize the \code{inputsTest} member.
#' if \code{dontNormTargets} is not set, then the targets are normalized in the same way. In classification problems,
#' normalizing the targets normally makes no sense. For regression, normalizing also the targets is usually a good idea. 
#' 
#' @title Function to normalize training and test set
#' @param x a list containing training and test data. Usually the output of \code{\link{splitForTrainingAndTest}}.
#' @param dontNormTargets should the target values also be normalized?
#' @param type type of the normalization. This parameter is passed to \code{\link{normalizeData}}. 
#' @return a named list with the same elements as \code{\link{splitForTrainingAndTest}}, but with normalized values.
#' The normalization parameters are appended to each member of the list as attributes, as in \code{\link{normalizeData}}.
#' @seealso \code{\link{splitForTrainingAndTest}}, \code{\link{normalizeData}}, \code{\link{denormalizeData}}, 
#' \code{\link{getNormParameters}}
#' @export
#' @examples
#' data(iris)
#' #shuffle the vector
#' iris <- iris[sample(1:nrow(iris),length(1:nrow(iris))),1:ncol(iris)]
#' 
#' irisValues <- iris[,1:4]
#' irisTargets <- decodeClassLabels(iris[,5])
#' 
#' iris <- splitForTrainingAndTest(irisValues, irisTargets, ratio=0.15)
#' normTrainingAndTestSet(iris)
normTrainingAndTestSet <- function(x, dontNormTargets=TRUE, type="norm") {
  
  inputsTrain <- normalizeData(x$inputsTrain, type=type)
  inputsTest <- normalizeData(x$inputsTest, type=attr(inputsTrain, "normParams"))
  
  if(dontNormTargets) {
    
    targetsTrain <- x$targetsTrain
    targetsTest <- x$targetsTest
    
  } else {
    
    targetsTrain <- normalizeData(x$targetsTrain, type=type)
    targetsTest <- normalizeData(x$targetsTest, type=attr(targetsTrain, "normParams"))
    
  }
  
  list(inputsTrain=inputsTrain, targetsTrain=targetsTrain, inputsTest=inputsTest, targetsTest=targetsTest)  
}



## Check the input for eventual problems.
checkInput <- function(x,y) {
  
  ok <- TRUE
  
  if(any(is.na(x))) {
    stop("missing values in 'x'")
    ok <- FALSE
  }
  if(any(is.na(y))) {
    stop("missing values in 'y'")
    ok <- FALSE
  }
  if(dim(x)[1L] != dim(y)[1L]) {
    stop("nrows of 'x' and 'y' must match")
    ok <- FALSE    
  }
  
  ok
}

#' Plot the iterative training and test error of the net of this \code{\link{rsnns}} object.
#'
#' Plots (if present) the class members \code{IterativeFitError} (as black line) and 
#' \code{IterativeTestError} (as red line).
#' 
#' @title Plot iterative errors of an rsnns object
# @param object the object to which to apply plotIterativeError
# @param ... additional function parameters
#' @export
plotIterativeError <- function(object, ...) UseMethod("plotIterativeError")

#' Plot the iterative training and test error of the net of this rsnns object.
#' 
#' @param object a rsnns object
#' @param ... parameters passed to \code{plot}
#' @export
# @S3method plotIterativeError rsnns
#' @method plotIterativeError rsnns
#' @rdname plotIterativeError
plotIterativeError.rsnns <- function(object, ...)
{
  if(!inherits(object, "rsnns")) stop("not a legitimate rsnns model")
  
  if(object$computeIterativeError) {
    plot(object$IterativeFitError, ylab="Weighted SSE", xlab="Iteration", type="l", ...)
    
    if(!is.null(object$IterativeTestError)) {
      
      testSetRatio <- nrow(as.matrix(object$fitted.values)) / nrow(as.matrix(object$fittedTestValues)) 
      
      lines(object$IterativeTestError * testSetRatio, col="red")    
    }    
  } else {
    print("Iterative error was not computed during training")
  }
  
}


#' This method decodes class labels from a numerical or levels vector to a
#' binary matrix, i.e.,  it converts the input vector to a binary matrix.
#' 
#' In the matrix, the value \code{valTrue} (e.g. 1) is present exactly in the
#' column given by the value in the input vector, and the value \code{valFalse}
#' (e.g. 0) in the other columns. The number of columns of the resulting matrix
#' depends on the number of unique labels found in the vector. E.g. the input
#' c(1, 3, 2, 3) will result in an output matrix with rows: 100 001 010 001
#' 
#' @title Decode class labels to a binary matrix
#' @references
#' Venables, W. N. and Ripley, B. D. (2002), 'Modern Applied Statistics with S', Springer-Verlag.
#' @param x class label vector
#' @param valTrue see Details paragraph
#' @param valFalse see Details paragraph
#' @return a matrix containing the decoded class labels
#' @export
#' @author The implementation is a slightly modified version of the function
#' \code{class.ind} from the \code{nnet} package of Brian Ripley.
#' @examples
#' decodeClassLabels(c(1,3,2,3))
#' decodeClassLabels(c("r","b","b","r", "g", "g"))
#' 
#' data(iris)
#' decodeClassLabels(iris[,5])
decodeClassLabels <- function(x, valTrue=1, valFalse=0)
{
  n <- length(x)
  x <- as.factor(x)
  res <- matrix(valFalse, n, length(levels(x)) )
  res[(1:n) + n*(unclass(x)-1)] <- valTrue
  dimnames(res) <- list(names(x), levels(x))
  
  res
}


#' Applies \code{analyzeClassification} row-wise to a matrix.
#' 
#' @title Encode a matrix of (decoded) class labels
#' @param x inputs
#' @param method see \code{analyzeClassification}
#' @param l idem
#' @param h idem
#' @return a numeric vector, each number represents a different class. A zero means
#' that no class was assigned to the pattern. 
#' @export
#' @seealso \code{\link{analyzeClassification}}
#' @examples 
#' data(iris)
#' labels <- decodeClassLabels(iris[,5])
#' encodeClassLabels(labels)
encodeClassLabels <- function(x, method="WTA", l=0.0, h=0.0) {
  apply(x, 1, function(y) analyzeClassification(y, method, l, h))
}

#' This function converts a vector (of class labels) to a numeric vector.
#' 
#' @title Convert a vector (of class labels) to a numeric vector
#' @param x inputs
#' @return the vector converted to a numeric vector
#' @export
#' @examples 
#' data(iris)
#' toNumericClassLabels(iris[,5])
toNumericClassLabels <- function(x) {
  if(is.numeric(x)) return(x)
  else return(as.numeric(x))
}

 
#' This function converts the continuous outputs to binary outputs that can be
#' used for classification. The two methods 402040, and winner-takes-all (WTA),
#' are implemented as described in the SNNS User Manual 4.2.
#' 
#' The following text is an edited citation from the SNNS User Manual 4.2 (pp
#' 269):
#' 
#' \describe{
#' \item{402040}{ 
#' A pattern is recognized as classified correctly, if (i) the output of exactly one output unit is >= h 
#' (ii) the teaching output of this unit is the maximum teaching output (> 0) of the pattern (iii) the output of all
#' other output units is <= l.
#' 
#' A pattern is recognized as classified incorrectly, if (i) and (iii) hold as above, but for (ii) holds that the teaching output is \emph{not}
#' the maximum teaching output of the pattern or there is no teaching output > 0.
#' 
#' A pattern is recognized as unclassified in all other cases.
#' 
#' The method derives its name from the commonly used default values l = 0.4, h = 0.6.
#' }
#' \item{WTA}{
#' A pattern is recognized as classified correctly, if (i) there is an output unit with the value greater than the output value of all other
#' output units (this output value is supposed to be a) (ii) a > h (iii) the teaching output of this unit is the maximum teaching output 
#' of the pattern (> 0) (iv) the output of all other units is < a - l.
#' 
#' A pattern is recognized as classified incorrectly, if (i), (ii), and (iv) hold as above, but for (iii) holds that the teaching output of this 
#' unit is \emph{not} the maximum teaching output of the pattern or there is no teaching output > 0.
#'  
#' A pattern is recognized as unclassified in all other cases. 
#' 
#' Commonly used default values for this method are: l = 0.0, h = 0.0.
#' }
#' } 
#' 
#' @title Converts continuous outputs to class labels
#' @param y inputs
#' @param method "WTA" or "402040"
#' @param l lower bound, e.g. in 402040: l=0.4
#' @param h upper bound, e.g. in 402040: h=0.6
#' @return the position of the winning unit (i.e., the winning class), or zero, if no classification was done.
#' @references 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen.
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' @export
#' @seealso \code{\link{encodeClassLabels}}
analyzeClassification <- function(y, method="WTA", l=0.0, h=0.0) {
  
  classes <- length(y)
  resClass <- 0
  
  if(method=="402040") {
    candClass <- which(y >= h)
    
    if(length(candClass) == 1) {
      if(max(y[-candClass]) <= l) {
        resClass <- candClass
      }
    }
    
  } else if(method=="WTA") {
    candClass <- which(y == max(y))
    
    if(length(candClass) == 1) {
      if(y[candClass] > h) {
        if(max(y[-candClass]) < (max(y) - l)) {
          resClass <- candClass
        }        
      }
    }
  } #else if(method=="ForcedWTA") {
    #resClass <- which(y==max(y))[1]
  #}
  
  resClass  
}

#' The confusion matrix shows how many times a pattern
#' with the real class x was classified as class y. A perfect method
#' should result in a diagonal matrix. All values not on the diagonal
#' are errors of the method.
#' 
#' If the class labels are not already encoded, they are encoded using \code{\link{encodeClassLabels}} 
#' (with default values). 
#' 
#' @title Computes a confusion matrix
#' @param targets the known, correct target values
#' @param predictions the corresponding predictions of a method for the targets
#' @return the confusion matrix
#' @export
confusionMatrix <- function(targets, predictions) {
  
  if(is.matrix(targets)) {
    if(ncol(targets)!=1) tr <- encodeClassLabels(targets)
  } else {
    tr <- toNumericClassLabels(targets)  
  }
  
  enc <- FALSE
  if(is.matrix(predictions)) {
    if(ncol(predictions)!=1) {
      pr <- encodeClassLabels(predictions)
      enc <- TRUE
    }
  } 
  
  if(!enc) pr <- toNumericClassLabels(predictions)
  
  cm <- table(targets=tr, predictions=pr)
  #rowSums(cm)
  #colSums(cm)
  cm
}


#' This function plots a receiver operating characteristic (ROC) curve. 
#' 
#' @title Plot a ROC curve
#' @param T predictions
#' @param D targets
#' @param ... parameters passed to plot
#' @references 
#' R news Volume 4/1, June 2004
#' @export
#' @author  Code is taken from R news Volume 4/1, June 2004.
plotROC <-function(T, D, ...){
  cutpoints<-c(-Inf, sort(unique(T)), Inf)
  sens<-sapply(cutpoints,
      function(c) sum(D[T>c])/sum(D))
  spec<-sapply(cutpoints,
      function(c) sum((1-D)[T<=c]/sum(1-D)))
  plot(1-spec, sens, type="l", ...)
}



#' The plot shows target values on the x-axis and fitted/predicted values on the y-axis. 
#' The optimal fit would yield a line through zero with gradient one.
#' This optimal line is shown in black color. A linear fit to the actual data
#' is shown in red color.
#' 
#' @title Plot a regression error plot
#' @param targets the target values
#' @param fits the values predicted/fitted by the model
#' @param ... parameters passed to \code{plot}
#' @export
plotRegressionError <- function(targets, fits, ...)
{
  
  plot(targets, fits, xlim=c(0,1), ylim=c(0,1), ...)
  
  linMod <- lm(fits ~ targets)
  abline(linMod, col="red")
  lines(c(0,1), c(0,1))
  
}
