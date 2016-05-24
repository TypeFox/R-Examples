# DESP/R/DESP_PEN_grad.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

DESP_PEN_grad <-
function(S,B,init,kappa,thresh,stepsize,tol){

  # read the sample size and the number of variables
  D = dim(S);
  p = D[2];               # p is the dimension

  SnB = diag(crossprod(S,B));

  f <- function(v) sum(-log(v) + SnB * v) + kappa * sum((B*t(B) > thresh) * (crossprod(diag(v),t(B)) - crossprod(t(B),diag(v)))^2 / (B^2 + t(B)^2), na.rm=TRUE)/2 

  df <- function(v) -1/v + SnB + kappa * 2*colSums((B*t(B) > thresh) * (crossprod(diag(v),t(B)) - crossprod(t(B),diag(v))) * (-B) / (B^2 + t(B)^2),na.rm=TRUE)
 
  x_prec = init
  x = x_prec+1;
  end = sum(x!=x_prec)!=0
  grad = df(x_prec);
  fx_prec = f(x_prec)

  stepsize=min(1,stepsize) # this line aims to prevent the stepsize from being to large, as a too large stepsize could make the gradient descent diverge
  i=1
  if (is.na(crossprod(grad))) stop("NA error in DESP_PEN_grad")
  while (i <= 5000 && end && crossprod(grad) >= tol) {

    # the multiplicative factors used for adaptative stepsize are those propose by Riedmiller and Braun in Rprop
    x <- x_prec - stepsize * grad/sqrt(crossprod(grad));
    x <- ifelse(x<1,1,x)
    end = sum(x!=x_prec)!=0
    fx = f(x)
    if(fx < fx_prec){
      stepsize = 1.2*stepsize
      grad = df(x);
      x_prec=x; 
      fx_prec = fx
    }else{
      stepsize = stepsize/2;
    }
    i=i+1
    if (is.na(crossprod(grad))) stop("NA error in DESP_PEN_grad (gradient possibly too large).")
  }

  return(x)
}
