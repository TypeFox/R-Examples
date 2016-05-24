# DESP/R/DESP_PML.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_PML <-
function(X,B,thresh,kappa,tol,Theta=NULL) {
  # estimation of the diagonal of the precision matrix by penalized likelihood minimization, when the true value of B is known or has already been estimated
  # the observations of the data matrix X are assumed to have zero mean
  # main function

  # read the sample size and the number of variables
  D = dim(X);
  n = D[1];               # n is the sample size
  p = D[2];               # p is the dimension

  hessian <-
  function(v,B,thresh,kappa) {
    p = dim(B)[2];            # p is the dimension
    H = matrix(0,p,p)

    for (i in 1:p){
      for (j in 1:p){
        if (B[j,i]*B[i,j]>thresh) {
          H[i,j] = kappa * 2 * (-B[i,j]) * B[j,i] / (B[i,j]^2 + B[j,i]^2)
    }}}

    for (j in 1:p){
      H[j,j] = 1/v[j]^2 
      for (j in 1:p){
        if (B[j,i]*B[i,j]>thresh) {
          H[j,j] = H[j,j] + kappa * 2 * B[i,j]^2 / (B[i,j]^2 + B[j,i]^2)
    }}}
    return(H)
  }

  # compute the sample cov matrix
  if(is.null(Theta))
    {
    S = crossprod(X)/n;
    }
  else
    {
    S = crossprod(X - Theta %*% MASS::ginv(B))/n;
    }
  
  # stepsize computation
  maxSV = max(svd(hessian(rep(n^(1/2),p),B,thresh,kappa))$d)
  
  Phi_inv = DESP_PEN_grad(S,B,rep(1,p),kappa,thresh,1/maxSV,tol);

  return(Phi_inv);
}
