# DESP/R/desp.cv.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

desp.cv <- function(X, v, lambda.range, gamma.range, settings=NULL){
  settings$outParCorrB <- TRUE
  n = nrow(X)
  p = ncol(X)

  minCrit <- Inf
  sn <- sample(n)
  subsets <- split(sn[1:(v*floor(n/v))], rep(1:v, each = floor(n/v)))
  for(lambda in lambda.range){
    for(gamma in gamma.range){
      crit <- 0
      for(i in 1:v){
        Xt <- X[setdiff(sn, subsets[[i]]),] # training
        Xv <- X[subsets[[i]],] # validation
        sSize <- length(subsets[[i]])
        pX <- desp(Xt, lambda, gamma, settings)
        crit <- crit + sum(sqrt(apply((Xv %*% pX$B - tcrossprod(rep(1,sSize), 
          crossprod(pX$B, pX$mu)))^2,1,sum)))/sSize
      }
      if(crit<minCrit){
        minCrit <- crit
        best <- c(lambda, gamma)
      }
    }
  }
  pX <- desp(X, best[1], best[2], settings)
  res <- list(Omega=pX$Omega, mu=pX$mu, Theta=pX$Theta, lambda=best[1], gamma=best[2])
  class(res) <- "desp.cv"
  return(res)
}
