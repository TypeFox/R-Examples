## Fri Mar 07 18:39:01 2014
## Original file Copyright © 2016 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################

######
######


OU <- function(N, ...)  UseMethod("OU")

OU.default <- function(N =100,M=1,x0=2,t0=0,T=1,Dt,mu=4,sigma=0.2,...)
               {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.numeric(sigma) || sigma <= 0) ) stop(" 'sigma' must be a numeric > 0 ")
    if (any(!is.numeric(mu) || mu <= 0) ) stop(" 'mu' must be a numeric > 0 ")
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    X <- HWV(N,M,x0,t0,T,Dt,mu,theta=0,sigma)
    return(X)
}
