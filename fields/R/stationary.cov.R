# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"stationary.cov" <- function(x1, x2=NULL, Covariance = "Exponential", Distance = "rdist", 
                             Dist.args = NULL, theta = 1, V = NULL, C = NA, marginal = FALSE, 
                             derivative = 0, distMat = NA, onlyUpper = FALSE, ...) {
  
  # get covariance function arguments from call
  cov.args <- list(...)
  # coerce x1 and x2 to matrices
  if (is.data.frame(x1)) 
    x1 <- as.matrix(x1)
  if (!is.matrix(x1)) 
    x1 <- matrix(c(x1), ncol = 1)
  
  if(is.null(x2))
    x2 = x1
  if (is.data.frame(x2)) 
    x2 <- as.matrix(x2)
  if (!is.matrix(x2)& !is.null(x2)) 
    x2 <- matrix(c(x2), ncol = 1)
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  #
  # separate out a single scalar transformation and a
  # more complicated scaling and rotation.
  # this is done partly to have the use of great circle distance make sense
  # by applying the scaling  _after_ finding the distance.
  #
  if (length(theta) > 1) {
    stop("theta as a vector matrix has been depreciated use the V argument")
  }
  #
  # following now treats V as a full matrix for scaling and rotation.
  #
  # try to catch incorrect conbination  of great circle distance and V
  if (Distance == "rdist.earth" & !is.null(V)) {
    stop("V not supported with great circle distance")
  }
  if (!is.null(V)) {
    if (theta != 1) {
      stop("can't specify both theta and V!")
    }
    x1 <- x1 %*% t(solve(V))
    x2 <- x2 %*% t(solve(V))
  }
  #
  # locations are now scaled and rotated correctly
  # now apply covariance function to pairwise distance matrix, or multiply
  # by C vector or just find marginal variance
  # this if block finds the cross covariance matrix
  if (is.na(C[1]) && !marginal) {
    #
    # if distMat is supplied, evaluate covariance for upper triangular part only
    #
    if(is.na(distMat[1])) {
      # distMat not supplied so must compute it along with covariance matrix
      # note overall scalling by theta (which is just theta under isotropic case)
      if(is.null(x2))
        distMat <- do.call(Distance, c(list(x1), Dist.args))
      else
        distMat <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
      
    }
    
    #
    # now convert distance matrix to covariance matrix
    #
    if(inherits(distMat, "dist")) {
      #
      # distMat is in compact form, so evaluate covariance over all distMat and convert to matrix form
      
      diagVal = do.call(Covariance, c(list(d=0), cov.args))
      
      if(onlyUpper)
        return(compactToMat(do.call(Covariance, c(list(d=distMat*(1/theta)), cov.args)), diagVal))
      else
        # if onlyUpper==FALSE, also set lower triangle of covariance matrix
        return(compactToMat(do.call(Covariance, c(list(d=distMat*(1/theta)), cov.args)), diagVal, lower.tri=TRUE))
    }
    else {
      # distMat is a full matrix
      return(do.call(Covariance, c(list(d = distMat/theta), cov.args)))
    }
  }
  # or multiply cross covariance by C
  # as coded below this is not particularly efficient of memory
  else if(!is.na(C[1])) {
    if(onlyUpper) {
      #the onlyUpper option is not compatible with the C option
      onlyUpper = FALSE
    }
    
    if(is.null(x2))
      bigD <- do.call(Distance, c(list(x1, x1), Dist.args))
    else
      bigD <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
    
    if (derivative == 0) {
      return(do.call(Covariance, c(list(d = bigD*(1/theta)), cov.args)) %*% C)
    }
    else {
      # find partial derivatives
      tempW <- do.call(Covariance, c(list(d = bigD*(1/theta)), 
                                     cov.args, derivative = derivative))
      # loop over dimensions and knock out each partial accumulate these in
      # in temp
      temp <- matrix(NA, ncol = d, nrow = n1)
      for (kd in 1:d) {
        # Be careful if the distance (tempD) is close to zero.
        # Note that the x1 and x2 are in transformed ( V inverse) scale
        sM <- ifelse(bigD == 0, 0, tempW * outer(x1[, kd], x2[, kd], "-")/(theta * bigD))
        # accumlate the new partial
        temp[, kd] <- sM %*% C
      }
      # transform back to original coordinates.
      if (!is.null(V)) {
        temp <- temp %*% t(solve(V))
      }
      return(temp)
    }
  }
  # or find marginal variance and return  a vector.
  else if (marginal) {
    sigma2 <- do.call(Covariance, c(list(d = 0), cov.args))
    return(rep(sigma2, nrow(x1)))
  }
  
  # should not get here based on sequence of conditional if statements above.
}
