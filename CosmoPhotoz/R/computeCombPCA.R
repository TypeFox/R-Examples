#  R package CosmoPhotoz file R/computeCombPCA.R
#  Copyright (C) 2014  COIN
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @title Combined PCA for training and test sample
#'
#' @description \code{computeCombPCA} computes combined PCA projections of 
#' the training and test samples.
#' 
#' @import pcaPP mvtnorm
#' @param x a matrix or a data.frame
#' @param y a matrix or a data.frame 
#' @param robust a boolean indicating if robust PCA should be used or not
#' @return PCA projections for each matrix 
#' @examples
#' 
#' #Multivariate data with outliers
#' library(mvtnorm)
#' x <- rbind(rmvnorm(100, rep(0, 6), diag(c(5, rep(1,5)))),
#'           rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
#' y <- rbind(rmvnorm(100, rep(0, 6), diag(c(5, rep(1,5)))),
#'           rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))         
#' #Here we calculate the principal components
#' pc <- computeCombPCA(x, y)
#'  
#' @usage computeCombPCA(x, y, robust)
#' 
#' @author Rafael S. de Souza, Alberto Krone-Martins
#' 
#' @keywords misc
#' @details The program is a simple alteration of PCAgrid() that computes a desired number of robust
#'  principal components using the grid search algorithm in the plane. 
#' @export
computeCombPCA <- function(x, y, robust=TRUE) {

  # First some basic error control
  if(!is.matrix(x)&!is.data.frame(x)) {
    stop("Error in computeCombPCA :: x is nor a matrix neither a data frame. The code expects a matrix or data frame.")
  }
  if(!is.matrix(y)&!is.data.frame(y)) {
    stop("Error in computeCombPCA :: y is nor a matrix neither a data frame. The code expects a matrix or data frame.")
  }
  
  # Now for the real work
  XY <- rbind(x, y)                                                # Create the combined matrix
  if(robust) {
    XYPCA <- PCAgrid(XY, k=ncol(XY)-1, scale="mad", method="mad")  # Calculate the robust PCA
    X.PCA <- as.data.frame(XYPCA$scores[1:nrow(x),])               # get the scores in a data frame
    Y.PCA <- as.data.frame(XYPCA$scores[(nrow(x)+1):nrow(XY),])    # get the scores in a data frame
   # PCvar <- which(cumsum((XYPCA$sdev)^2) / sum(XYPCA$sdev^2) >= npcvar)[1]
  } else {
    XYPCA <- prcomp(XY, scale=TRUE)                                # non robust PCA
    XYPCA$x <- XYPCA$x[, 1:ncol(XY)-1]                             # get the same columns as the robust
    nn <- c()                                                      # create the name of the columns
    for(i in 1:ncol(XYPCA$x)) { 
      nn <- c(nn, paste("Comp.",i,sep=""))
    }
    X.PCA <- as.data.frame(XYPCA$x[1:nrow(x),])                    # split and get the scores in a data frame
    Y.PCA <- as.data.frame(XYPCA$x[(nrow(x)+1):nrow(XY),])         # split and get the scores in a data frame 
    names(X.PCA) <- nn                                             # set the column names as the robust
    names(Y.PCA) <- nn
  }
  
  # That's all folks!
  return(list(x=X.PCA, y=Y.PCA, PCsum=summary(XYPCA)))
}

