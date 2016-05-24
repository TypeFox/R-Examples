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


BB <- function(N, ...)  UseMethod("BB")

BB.default <- function(N =100,M=1,x0=0,y=1,t0=0,T=1,Dt,...)
             {
    if (any(!is.numeric(x0) || !is.numeric(y)) ) stop("'x0' and 'y' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(t0 < 0 || T < 0 || T <= t0) ) 
        stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N
    bb <- function(){
         w = c(0,cumsum(rnorm(N,mean=0,sd=sqrt(Dt))))
         X <- x0 + w - (t-t0)/(T-t0) * (w[N+1]-y +x0)
         X
         }
    res <- data.frame(sapply(1:M,function(i) bb()))
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = t0, deltat = Dt)
    return(X)
}



