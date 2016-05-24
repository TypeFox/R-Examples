#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/




networkInferenceGlassoBIC <- function(dataNet, nb.rho = 100){
  

  ## perform network inference using the glasso 
  ## from data and select best model using bic
  ## input dataNet (matrix or data frame): data used to perform network inference
  ## input nb.rho (integer): number of regularization parameters rho tested
  ## output A (matrix): selected adjacency matrix based on BIC
  ## output Theta (matrix): selected theta matrix based on BIC
  ## output Sigma (matrix): selected Sigma matrix based on BIC
  ## outpout penaltieslist (list): list of regularization parameters
  ## outpout pathA (list): list of adjacency matrix for each regularization parameter
  ## outpout pathTheta (list): list of sigma theta for each regularization parameter
  ## outpout pathSigma (list): list of sigma matrix for each regularization parameter
  ## outpout pathBIC (list): list of BIC value for each regularization parameter
 
  ## require 
  ## library(glasso)
  
  if(is.matrix(dataNet) == FALSE & is.data.frame(dataNet) == FALSE) 
    stop(paste(sQuote("dataNet"), "must be a matrix"))
  
  ## data scaling
  dataNet <- scale(dataNet)
  ## compute an appropriate list of regularization parameters
  ## to use as an input of the glasso function
  tasks = factor(rep(1,nrow(dataNet)))
  penalty.min <- 1e-02
  penalty.max <- max(sapply(sapply(by(dataNet, tasks, var), abs), max))
  penalty.min <- max(penalty.min, .Machine$double.eps)
  penaltieslist <- seq(penalty.max, penalty.min, length = nb.rho)
  ## compute sample covariance matrix for input in the glassopath function
  forglasso <- cov(dataNet)
  resGlasso <- glassopath(forglasso, rholist=penaltieslist, penalize.diagonal=FALSE, trace=0)
  ## collect results from the glasso output
  pathTheta <- lapply(seq(dim(resGlasso$wi)[3]), function(x) resGlasso$wi[ , , x])
  pathSigma <- lapply(seq(dim(resGlasso$w)[3]), function(x) resGlasso$w[ , , x])
  pathA <-pathTheta
  pathA <- lapply(pathA,function(mat){ mat[mat!=0] <- 1; return(mat)} )
  ## select a network based on the BIC criterion
  Sbic <- var(dataNet, na.rm = TRUE)
  getBIC <- function(dumTheta){
    loglikfun <- (dim(dataNet)[1]/2) * (log(det(dumTheta)) - sum(diag(dumTheta %*% Sbic)))
    dffun <- sum(abs(dumTheta) > 0) - sum(abs(diag(dumTheta)) > 0)
    BIC <- loglikfun - (log(dim(dataNet)[1])/2)*dffun/2
    return(BIC)
  }
  pathBIC <- unlist(lapply(pathTheta, getBIC ))
  return(list(A = pathA[[which.max(pathBIC)]], Sigma = pathSigma[[which.max(pathBIC)]], penaltieslist = penaltieslist, Theta=pathTheta[[which.max(pathBIC)]], resGlasso = resGlasso, pathA = pathA, pathSigma = pathSigma, pathTheta = pathTheta, pathBIC = pathBIC))
}
