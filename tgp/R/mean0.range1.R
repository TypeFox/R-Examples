#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## mean0.range1:
##
## translate the input columns (X.m) to each (independently)
## have a mean of zero and a range of one -- as used by Chipman
## et al.  Also save the necessary mean and range information
## so that the transformation can be undone later

"mean0.range1" <-
function(X.m) {

  ## checks and coersion into a matrix
  if(is.null(X.m)) return(NULL)
  else if(is.null(dim(X.m))) X <- matrix(X.m, ncol=1)
  else X <- X.m

  ## initialize the information necesary for undoing
  undo <- list()
  undo$min <- rep(0, ncol(X))
  undo$max <- rep(0, ncol(X))
  undo$amean <- rep(0, ncol(X))

  ## make the transformation in each dimension
  for(i in 1:ncol(X)) {
    undo$min[i] <- min(X[,i])
    undo$max[i] <- max(X[,i])
    X[,i] <- X[,i] / (max(X[,i]) - min(X[,i]))
    undo$amean[i] <- mean(X[,i])
    X[,i] <- X[,i] - mean(X[,i])
  }

  ## convert input vectors back into vectors
  if(is.null(dim(X.m))) X.m <- as.vector(X)
  else X.m <- X

  ## return both the transformed data and the info to undo
  return(list(X=X,undo=undo))
}


## undo.mean0.range1:
##
## using the info saved by mean0.range1, undo the transformation
## on X -- usually the undo is performed on new data that is curently
## on the scale of the transformed X, but should be reported on the
## scale of the original (unransformed) X

"undo.mean0.range1" <- 
function(X.m, undo, nomean=FALSE, s2=FALSE)
{
  ## checks and coerse into a matrix
  if(is.null(X.m)) return(NULL)
  else if(is.null(dim(X.m))) X <- matrix(X.m, ncol=1)
  else X <- X.m

  ## undo in each column of X
  for(i in 1:(dim(X)[2])) {
    if(!nomean) X[,i] <- X[,i] + undo$amean[i]
    if(s2) X[,i] <- X[,i]*(undo$max[i] - undo$min[i])^2
    else X[,i] <- X[,i]*(undo$max[i] - undo$min[i])
  }

  ## convert input vectors back into vectors
  if(is.null(dim(X.m))) X.m <- as.vector(X)
  else X.m <- X

  ## return the undone transformation
  return(X.m)
}
