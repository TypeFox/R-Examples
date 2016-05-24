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


## check.matrix:
##
## check/enfore that the X matrix has the proper dimensions
## (and the same number or rows as length(Z)) removing invalid rows
## (of Z too), i.e., NA, NaN, Inf

"check.matrix" <- 
function(X, Z=NULL)
{
  ## format X
  if(is.null(X)) return(NULL)
  n <- nrow(X)
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n) }
  X <- as.matrix(X)
  
  ## if a Z is provided to go along with X
  if(!is.null(Z)) {
    
    ## format Z
    ## Z <- as.vector(matrix(Z, ncol=1)[,1])
    Z <- as.vector(as.matrix(Z))
    if(length(Z) != n) stop("mismatched row dimension in X and Z")
    
    ## calculate locations of NAs NaNs and Infs in Z
    nna <- (1:n)[!is.na(Z) == 1]
    nnan <- (1:n)[!is.nan(Z) == 1]
    ninf <- (1:n)[!is.infinite(Z) == 1]
    if(length(nna) < n) warning(paste(n-length(nna), "NAs removed from input vector"))
    if(length(nnan) < n) warning(paste(n-length(nnan), "NaNs removed from input vector"))
    if(length(ninf) < n) warning(paste(n-length(ninf), "Infs removed from input vector"))
    
    neitherZ <- intersect(nna, intersect(nnan, ninf))
  } else neitherZ <- (1:n)
  
  ## calculate row locations of NAs NaNs and Infs in X
  nna <- (1:n)[apply(!is.na(X), 1, prod) == 1]
  nnan <- (1:n)[apply(!is.nan(X), 1, prod) == 1]
  ninf <- (1:n)[apply(!is.infinite(X), 1, prod) == 1]
  if(length(nna) < n) warning(paste(n-length(nna), "NAs removed from input matrix"))
  if(length(nnan) < n) warning(paste(n-length(nnan), "NaNs removed from input matrix"))
  if(length(ninf) < n) warning(paste(n-length(ninf), "Infs removed from input matrix"))
  neitherX <- intersect(nna, intersect(nnan, ninf))

  ## oops, no data:
  if(length(neitherX) == 0)
    stop("no valid (non-NA NaN or Inf) data found")

  ## combine good X and Z rows
  neither <- intersect(neitherZ, neitherX)
  X <- matrix(X[neither,], nrow=length(neither))
  Z <- Z[neither]

  return(list(X=X, Z=Z))
}


## famify.X
##
## change an X matrix into a data frame with the names specified
## used by tgp.postprocess to convert a matrix enforced by check.matrix
## back into the data frame it started as

"framify.X" <-
function(X, Xnames, d)
{
  X <- data.frame(t(matrix(X, nrow=d)))
  if(is.null(Xnames)) {
    nms <- c();
    for(i in 1:d) { nms <- c(nms, paste("x", i, sep="")) }
    names(X) <- nms
  } else { names(X) <- Xnames }
  return(X)
}

