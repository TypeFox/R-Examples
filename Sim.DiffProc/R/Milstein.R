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


#####
##### Milstein1D
 
.Milstein1D <- function(N =100,M=1,x0=1,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
                       {
    DSx  <- D(diffusion,"x")  
    if (type=="ito"){A    <- function(t,x)  eval(drift)}else{
    A  <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(DSx)}
    S  <- function(t,x)  eval(diffusion)
    Sx <- function(t,x)  eval(DSx)
    x0 <- rep(x0,M)[1:M]	
    X <- .Call("Milstein1d", x0, t0, Dt, as.integer(N), as.integer(M), A, S, Sx, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      

#####
##### Milstein2D

.Milstein2D <- function(N =100,M=1,x0=2,y0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
                       {					     
    DSx  <- D(diffx,"x")
    DSy  <- D(diffy,"y")  
    if (type=="ito"){Ax    <- function(t,x,y) eval(driftx)
                     Ay    <- function(t,x,y) eval(drifty)}else{
    Ax  <- function(t,x,y) eval(driftx) - 0.5 * eval(diffx) * eval(DSx)
    Ay  <- function(t,x,y) eval(drifty) - 0.5 * eval(diffy) * eval(DSy)}
    Sx  <- function(t,x,y) eval(diffx)
    Sy  <- function(t,x,y) eval(diffy) 
    dSx <- function(t,x,y) eval(DSx)
    dSy <- function(t,x,y) eval(DSy) 
    x0 <- rep(x0,M)[1:M]
    y0 <- rep(y0,M)[1:M]	
    Val <- .Call("Milstein2d", x0, y0, t0, Dt, as.integer(N), as.integer(M), Ax, Sx, dSx, Ay, Sy, dSy, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(Val[,1:M], start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Val[,(M+1):(2*M)], start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}


#####
##### Milstein3D

.Milstein3D <- function(N =100,M=1,x0=2,y0=1,z0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    DSx  <- D(diffx,"x")
    DSy  <- D(diffy,"y")
    DSz  <- D(diffz,"z")  
    if (type=="ito"){
        Ax    <- function(t,x,y,z) eval(driftx)
        Ay    <- function(t,x,y,z) eval(drifty)
        Az    <- function(t,x,y,z) eval(driftz)}else{
    Ax <- function(t,x,y,z) eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    Ay <- function(t,x,y,z) eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    Az <- function(t,x,y,z) eval(driftz) - 0.5 * eval(diffz) * eval(D(diffz,"z"))}
    Sx  <- function(t,x,y,z) eval(diffx)
    Sy  <- function(t,x,y,z) eval(diffy)
    Sz  <- function(t,x,y,z) eval(diffz) 
    dSx <- function(t,x,y,z) eval(DSx)
    dSy <- function(t,x,y,z) eval(DSy)
    dSz <- function(t,x,y,z) eval(DSz)
    x0 <- rep(x0,M)[1:M]
    y0 <- rep(y0,M)[1:M]
    z0 <- rep(z0,M)[1:M]
    Val <- .Call("Milstein3d", x0, y0, z0, t0, Dt, as.integer(N), as.integer(M), Ax, Sx, dSx, Ay, Sy , dSy, Az, Sz, dSz, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(Val[,1:M], start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Val[,(M+1):(2*M)], start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Val[,(2*M+1):(3*M)], start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 

