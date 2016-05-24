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


## friedman.1.data:
##
## generate a random sample of size n from Friedman's 10-d
## first data set used to validate the MARS method -- the
## response depends linearly and non-linearly on the first
## five inputs only

"friedman.1.data" <-
function(n=100)
{
  X <- matrix(runif(n*10), nrow=n)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3]-0.5)^2 + 10*X[,4] + 5*X[,5]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X,Y,Ytrue))
}


## fried.bool:
##
## generate a random sample of size n from a boolean segmented
## version of Friedman's 10-d first data set used to validate the
## MARS method -- the response depends linearly and non-linearly
## on the first five inputs only, but now which part of the function
## is on depends on an indicator 1:4

"fried.bool" <-
function(n=100)
{
  ## a function that is a sum of parts
  f1 <- function(X) { 10*sin(pi*X[,1]*X[,2]) }
  f2 <- function(X) { 20*(X[,3]-0.5)^2 }
  f3 <- function(X) { 10*X[,4] + 5*X[,5] }
  f4 <- function(X) {
    10*sin(pi*X[,5]*X[,4]) + 20*(X[,3]-0.5)^2 + 10*X[,2] + 5*X[,1]
  }
  fs <- list(f1, f2, f3, f4)

  ## boolean codings of 1:4
  BoolI <- rbind(c(0,0,0), c(0,0,1), c(0,1,0), c(1,0,0))

  ## sample n indicators in 1:4 and record their boolean coding
  I <- sample(c(1,2,3,4), n, replace=TRUE)
  Imat <-matrix(BoolI[I,], nrow=n)

  ## n random inputs U(0,1) in 10 dimensions
  X <- matrix(runif(n*10), nrow=n)

  ## allocate space for the true response
  Ytrue <- rep(0, n)

  ## calculate responses for  each of the four groups
  for(i in 1:4) {
    indx <- I == i
    if(sum(indx) == 0) next;
    indx <- (1:n)[indx]
    XX <- matrix(X[indx,], ncol=10)
    Ytrue[indx] <- fs[[i]](XX)
  }

  ## add some noise
  Y <- Ytrue + rnorm(n, 0, 1)

  ## return the inputs, bookean coding and outputs
  return(data.frame(X=X, I=Imat, Y, Ytrue))
}
