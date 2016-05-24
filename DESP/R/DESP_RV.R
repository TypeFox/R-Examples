# DESP/R/DESP_RV.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

DESP_RV <-
function(X,B,Theta=NULL) {
  # estimation of the diagonal of the precision matrix by residual variance when the true value of B is known or has already been estimated
  # the observations of the data matrix X are assumed to have zero mean
  
  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size

  normE2 <- function(x){
    # squared Euclidean norm of a vector
    sum(x^2)
  }

  if(is.null(Theta))
    {
    Phi <- apply(tcrossprod(X,t(B)),2,normE2)/n;
    }
  else
    {
    Phi <- apply(tcrossprod(X,t(B))-Theta,2,normE2)/n;
    }

  return(1/Phi);
}

DESP_AD <-
function(X,B,Theta=NULL) {
  # estimation of the diagonal of the precision matrix by average absolute deviation around the mean when the true value of B is known or has already been estimated
  # the observations of the data matrix X are assumed to have zero mean
  
  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size

  norm1 <- function(x){
    # l1 norm of a vector
    sum(abs(x))
  }

  if(is.null(Theta))
    {
    Phi <- apply(tcrossprod(X,t(B)),2,norm1)^2 * pi/2 /n^2;
    }
  else
    {
    Phi <- apply(tcrossprod(X,t(B))-Theta,2,norm1)^2 * pi/2 /n^2;
    }

  return(1/Phi);
}

