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
##### STS1D

.STS1D <- function(N =100,M=1,x0=1,t0=0,T=1,Dt,drift,diffusion,
                          type=c("ito","str"),...)
                       {
	DSx  <- D(diffusion,"x")  
    DSxx <- D(DSx,"x")
    if (type=="ito"){
    A    <- function(t,x)  eval(drift)
    Ax   <- function(t,x)  eval(D(drift,"x"))
    Axx  <- function(t,x)  eval(D(D(drift,"x"),"x"))
    }else{
    A    <- function(t,x)  eval(drift) - 0.5 * eval(diffusion) * eval(D(diffusion,"x"))
    Ax   <- function(t,x)  eval(D(drift,"x")) - 0.5 * (eval(D(diffusion,"x")) * eval(D(diffusion,"x"))+ eval(diffusion) * eval(D(D(diffusion,"x"),"x")))
    Axx  <- function(t,x)  eval(D(D(drift,"x"),"x")) - 0.5 * ( eval(D(D(diffusion,"x"),"x")) * eval(D(diffusion,"x"))+ eval(D(diffusion,"x")) * eval(D(D(diffusion,"x"),"x"))+
                           eval(D(diffusion,"x")) * eval(D(D(diffusion,"x"),"x")) + eval(diffusion) * eval(D(D(D(diffusion,"x"),"x"),"x")) )
                  }
    S    <- function(t,x)  eval(diffusion)
    Sx   <- function(t,x)  eval(DSx)
    Sxx  <- function(t,x)  eval(DSxx)
    x0 <- rep(x0,M)[1:M]	
	Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), 2, 2)
    RN <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    Z <- RN[, 1]
    U <- RN[, 2]
    X <- .Call("Sts1d", x0, t0, Dt, as.integer(N), as.integer(M), A, Ax, Axx, S, Sx, Sxx, Z, U, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- "X"
    name <- if(M > 1) paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=name)
    return(list(X=X))
}   
      
             

#####
##### STS2D

.STS2D <- function(N =100,M=1,x0=2,y0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                          type=c("ito","str"),...)
                       {
    DSx  <- D(diffx,"x")  
    DSxx <- D(DSx,"x")
    DSy  <- D(diffy,"y")  
    DSyy <- D(DSy,"y")					   
    if (type=="ito"){
    Ax    <- function(t,x,y)  eval(driftx)
    dAx   <- function(t,x,y)  eval(D(driftx,"x"))
    dAxx  <- function(t,x,y)  eval(D(D(driftx,"x"),"x"))
    Ay    <- function(t,x,y)  eval(drifty)
    dAy   <- function(t,x,y)  eval(D(drifty,"y"))
    dAyy  <- function(t,x,y)  eval(D(D(drifty,"y"),"y"))
    }else{
    Ax    <- function(t,x,y)  eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    dAx   <- function(t,x,y)  eval(D(driftx,"x")) - 0.5 * (eval(D(diffx,"x")) * eval(D(diffx,"x"))+ eval(diffx) * eval(D(D(diffx,"x"),"x")))
    dAxx  <- function(t,x,y)  eval(D(D(driftx,"x"),"x")) - 0.5 * ( eval(D(D(diffx,"x"),"x")) * eval(D(diffx,"x"))+ eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x"))+
                              eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x")) + eval(diffx) * eval(D(D(D(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y)  eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    dAy   <- function(t,x,y)  eval(D(drifty,"y")) - 0.5 * (eval(D(diffy,"y")) * eval(D(diffy,"y"))+ eval(diffy) * eval(D(D(diffy,"y"),"y")))
    dAyy  <- function(t,x,y)  eval(D(D(drifty,"y"),"y")) - 0.5 * ( eval(D(D(diffy,"y"),"y")) * eval(D(diffy,"y"))+ eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y"))+
                              eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y")) + eval(diffy) * eval(D(D(D(diffy,"y"),"y"),"y")) )
                  }
    Sx    <- function(t,x,y)  eval(diffx)
    dSx   <- function(t,x,y)  eval(DSx)
    dSxx  <- function(t,x,y)  eval(DSxx)
    Sy    <- function(t,x,y)  eval(diffy)
    dSy   <- function(t,x,y)  eval(DSy)
    dSyy  <- function(t,x,y)  eval(DSyy)
    x0 <- rep(x0,M)[1:M]
    y0 <- rep(y0,M)[1:M]
    Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), 2, 2)
    RN1 <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    RN2 <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    Z1 <- RN1[, 1]
    U1 <- RN1[, 2]	
    Z2 <- RN2[, 1]
    U2 <- RN2[, 2]	
    Val <- .Call("Sts2d", x0, y0, t0, Dt, as.integer(N), as.integer(M), Ax, dAx, dAxx, Ay, dAy, dAyy, Sx, dSx, dSxx, Sy, dSy,
                  dSyy, Z1, U1, Z2, U2, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- c("X","Y")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""))
    X <- ts(Val[,1:M], start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Val[,(M+1):(2*M)], start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    return(list(X=X,Y=Y))
}

#####
##### STS3D

.STS3D <- function(N =100,M=1,x0=2,y0=1,z0=1,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                     driftz,diffz,type=c("ito","str"),...)
                       {
    if (type=="ito"){
    Ax    <- function(t,x,y,z)  eval(driftx)
    dAx   <- function(t,x,y,z)  eval(D(driftx,"x"))
    dAxx  <- function(t,x,y,z)  eval(D(D(driftx,"x"),"x"))
    Ay    <- function(t,x,y,z)  eval(drifty)
    dAy   <- function(t,x,y,z)  eval(D(drifty,"y"))
    dAyy  <- function(t,x,y,z)  eval(D(D(drifty,"y"),"y"))
    Az    <- function(t,x,y,z)  eval(driftz)
    dAz   <- function(t,x,y,z)  eval(D(driftz,"z"))
    dAzz  <- function(t,x,y,z)  eval(D(D(driftz,"z"),"z"))}else{
    Ax    <- function(t,x,y,z)  eval(driftx) - 0.5 * eval(diffx) * eval(D(diffx,"x"))
    dAx   <- function(t,x,y,z)  eval(D(driftx,"x")) - 0.5 * (eval(D(diffx,"x")) * eval(D(diffx,"x"))+ eval(diffx) * eval(D(D(diffx,"x"),"x")))
    dAxx  <- function(t,x,y,z)  eval(D(D(driftx,"x"),"x")) - 0.5 * ( eval(D(D(diffx,"x"),"x")) * eval(D(diffx,"x"))+ eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x"))+
                                eval(D(diffx,"x")) * eval(D(D(diffx,"x"),"x")) + eval(diffx) * eval(D(D(D(diffx,"x"),"x"),"x")) )
    Ay    <- function(t,x,y,z)  eval(drifty) - 0.5 * eval(diffy) * eval(D(diffy,"y"))
    dAy   <- function(t,x,y,z)  eval(D(drifty,"y")) - 0.5 * (eval(D(diffy,"y")) * eval(D(diffy,"y"))+ eval(diffy) * eval(D(D(diffy,"y"),"y")))
    dAyy  <- function(t,x,y,z)  eval(D(D(drifty,"y"),"y")) - 0.5 * ( eval(D(D(diffy,"y"),"y")) * eval(D(diffy,"y"))+ eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y"))+
                                eval(D(diffy,"y")) * eval(D(D(diffy,"y"),"y")) + eval(diffy) * eval(D(D(D(diffy,"y"),"y"),"y")) )
    Az    <- function(t,x,y,z)  eval(driftz) - 0.5 * eval(diffz) * eval(D(diffz,"z"))
    dAz   <- function(t,x,y,z)  eval(D(driftz,"z")) - 0.5 * (eval(D(diffz,"z")) * eval(D(diffz,"z"))+ eval(diffz) * eval(D(D(diffz,"z"),"z")))
    dAzz  <- function(t,x,y,z)  eval(D(D(driftz,"z"),"z")) - 0.5 * ( eval(D(D(diffz,"z"),"z")) * eval(D(diffz,"z"))+ eval(D(diffz,"z")) * eval(D(D(diffz,"z"),"z"))+
                                eval(D(diffz,"z")) * eval(D(D(diffz,"z"),"z")) + eval(diffz) * eval(D(D(D(diffz,"z"),"z"),"z")) )
                  }
    DSx  <- D(diffx,"x")  
    DSxx <- D(DSx,"x")
    Sx    <- function(t,x,y,z)  eval(diffx)
    dSx   <- function(t,x,y,z)  eval(DSx)
    dSxx  <- function(t,x,y,z)  eval(DSxx)
    DSy  <- D(diffy,"y")  
    DSyy <- D(DSy,"y")
    Sy    <- function(t,x,y,z)  eval(diffy)
    dSy   <- function(t,x,y,z)  eval(DSy)
    dSyy  <- function(t,x,y,z)  eval(DSyy)
    DSz  <- D(diffz,"z")  
    DSzz <- D(DSz,"z")
    Sz    <- function(t,x,y,z)  eval(diffz)
    dSz   <- function(t,x,y,z)  eval(DSz)
    dSzz  <- function(t,x,y,z)  eval(DSzz)
    x0 <- rep(x0,M)[1:M]
    y0 <- rep(y0,M)[1:M]
    z0 <- rep(z0,M)[1:M]
    Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, (1/3)*Dt^3), 2, 2)
    RN1 <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    RN2 <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    RN3 <- .rMnorm(N*M, Mu= c(0, 0), Sigma)
    Z1 <- RN1[, 1]
    U1 <- RN1[, 2]	
    Z2 <- RN1[, 1]
    U2 <- RN1[, 2]	
    Z3 <- RN1[, 1]
    U3 <- RN1[, 2]	
    Val <- .Call("Sts3d", x0, y0, z0, t0, Dt, as.integer(N), as.integer(M), Ax,dAx, dAxx, Ay,dAy, dAyy, Az,dAz, dAzz,
                 Sx, dSx, dSxx, Sy, dSy, dSyy, Sz, dSz, dSzz, Z1, U1, Z2, U2, Z3, U3, .GlobalEnv, PACKAGE="Sim.DiffProc")
    name <- c("X","Y","Z")
    name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
    X <- ts(Val[,1:M], start = t0, deltat = Dt, names=name[1:M])
    Y <- ts(Val[,(M+1):(2*M)], start = t0, deltat = Dt, names=name[(M+1):(2*M)])
    Z <- ts(Val[,(2*M+1):(3*M)], start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
    return(list(X=X,Y=Y,Z=Z))
} 



