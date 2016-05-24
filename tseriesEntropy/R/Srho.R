
  ## tseriesEntropy
  ## Entropy based tests of serial dependence and nonlinearity
  #
  #  The authors of this software is
  #
  #  Simone Giannerini, Copyright (c) 2009.
  #
  #  Permission to use, copy, modify, and distribute this software for any
  #  purpose without fee is hereby granted, provided that this entire notice
  #  is included in all copies of any software which is or includes a copy
  #  or modification of this software and in all copies of the supporting
  #  documentation for such software.
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

## ***************************************************************************************************
setClass("Srho",
         representation(.Data="numeric",lags="integer",stationary="logical",data.type="character", notes = "character"))


## ***************************************************************************************************
setMethod("plot" , signature(x = "Srho",y = "missing"),
     function(x, y, type = "s", xlab = "lag", ylab = "S", ylim=c(0,max(x@.Data)), main = NULL,col=1,mai=c(.85,.75,.1,.1),lwd=1.5, grid=TRUE, ...){
            par(mai= mai);
            plot(x@lags,x@.Data,type=type,col=col,xlab=xlab,ylab=ylab,ylim=ylim,lwd=lwd,...);
            if(grid){grid()}
      }
)
## ***************************************************************************************************
setMethod ("show" , "Srho",
    function(object){
    out <- object@.Data
    names(out) <- object@lags;
    n <- length(out);
    lag.max <- object@lags[n]
    cat (" Srho computed on", lag.max, "lags \n")
    cat (" -------------------------------------------------------------------------- \n")
    print(out)
    cat (" -------------------------------------------------------------------------- \n")
    cat (" Data type          :" , object@data.type , "\n")
    cat (" Stationary version :" , object@stationary , "\n")
    if(length(object@notes)>0){
    cat (" Additional notes   :" ,   object@notes, "\n");
    }
    cat ("\n")
    }
)

## ***************************************************************************************************

Srho <- function(x,y,lag.max,stationary=TRUE, plot=TRUE,version=c("FORTRAN","R"), nor=FALSE){
    version <- match.arg(version);

    if(!(is.integer(x)||is.factor(x)||is.character(x))) stop('input series must be either integer or categorical valued')
    n <- length(x);
    if (missing(lag.max)) lag.max = round(n/4);
    if((lag.max >= n)||(lag.max < 2)) stop('incorrect number of lags')


    if (missing(y)){
        if(version =='FORTRAN') {
            return(Srho.uni.F(x,lag.max=lag.max,stationary=stationary, plot=plot, nor=nor))
        }else if(version=='R'){
            return(Srho.uni.R(x,lag.max=lag.max,stationary=stationary, plot=plot, nor=nor))
        }else {stop("A version must be specified: either FORTRAN or R")}
    }else{
        if(!(is.integer(y)||is.factor(y))) stop('input series must be either integer or categorical valued')
        if(version =='FORTRAN'){
            return(Srho.biv.F(x,y,lag.max=lag.max,stationary=stationary, plot=plot, nor=nor))
        }else if(version=='R'){
            return(Srho.biv.R(x,y,lag.max=lag.max,stationary=stationary, plot=plot, nor=nor))
        }else {stop("A version must be specified: either FORTRAN or R")}
    }
}
## ***************************************************************************************************

Srho.uni.F <- function(x,lag.max,stationary=TRUE, plot=FALSE, nor=FALSE){

    ## An implementation of the entropy-based dependence measure S[rho]
    ## proposed by Granger et al. (2004) Journal of Time Series Analysis.
    ## Deals with INTEGER or CATEGORICAL time series, for continuous processes non-parametric
    ## density estimation is required.
    ## Probabilities are estimated through relative frequencies

    ## PARAMETERS:
    ##   x:          Integer or Categorical series
    ##   lag.max:       number of lags  -- default: (round(n/4)
    ##   stationary :  logical - if TRUE computes the version assuming stationarity (default)

    ## Simone Giannerini  2007

    ## FORTRAN VERSION OF Srho.R

    if(stationary==TRUE) {
        Srho <- .Fortran("ssuni2",as.integer(x),as.integer(length(x)),as.integer(lag.max),S=as.double(rep(0,lag.max)), as.integer(nor),PACKAGE="tseriesEntropy")$S
    }
    else {
        Srho <- .Fortran("ssuni",as.integer(x),as.integer(length(x)),as.integer(lag.max),S=as.double(rep(0,lag.max)), as.integer(nor),PACKAGE="tseriesEntropy")$S
    }
    out <- new("Srho")
    out@.Data      <- Srho
    out@lags       <- 1:lag.max
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}
 ## ***************************************************************************************************

Srho.uni.R <- function(x,lag.max,stationary=TRUE, plot=FALSE, nor=FALSE){

    ## An implementation of the entropy-based dependence measure S[rho]
    ## proposed by Granger et al. (2004) Journal of Time Series Analysis.
    ## Deals with INTEGER or CATEGORICAL time series, for continuous processes non-parametric
    ## density estimation is required.
    ## Probabilities are estimated through relative frequencies

    ## PARAMETERS:
    ##   x          :  Integer or Categorical series
    ##   lag.max    :  number of lags  -- default: (round(n/4)
    ##   stationary :  logical - if TRUE computes the version assuming stationarity (default)

    ## Simone Giannerini  2007

    ## R VERSION OF Srho.F
    n <- length(x);
    Srho    <- double(lag.max)
    if(stationary==TRUE) {
        md   <- prop.table(table(x))
        i    <- dim(md);
        for(k in 1:lag.max) {
            x1  <- x[1:(n-k)];
            x2  <- x[(k+1):n];
            jd <- matrix(0,i,i)
            dum  <- prop.table(table(x1,x2));
            dimnames(jd) <- list(names(md),names(md))
            jd[dimnames(dum)$x1,dimnames(dum)$x2] <- dum
            S <- 0;
            for(k1 in 1:i) {
                for(k2 in 1:i) {
                    S <- S + (sqrt(jd[k1,k2])-sqrt(md[k1]*md[k2]))^2;
                }
            }
            Srho[k] <- S;
        }
        if(nor){
            smax <- 1 - sum(md^(3/2))
            Srho <- Srho/smax
        }
    }
    else {
        for(k in 1:lag.max) {
            x1  <- x[1:(n-k)];
            x2  <- x[(k+1):n];
            md1 <- table(x1);
            md2 <- table(x2);
            jd  <- table(x1,x2);
            i1  <- length(md1);
            i2  <- length(md2);
            S <- 0;
            for(k1 in 1:i1) {
                for(k2 in 1:i2) {
                    S <- S + (sqrt(jd[k1,k2]*(n-k))-sqrt(md1[k1]*md2[k2]))^2;
                }
            }
            Srho[k] <- (S/(n-k)^2);
            if(nor){
                smax <- max(1 - sum((md1/(n-k))^(3/2)),1 - sum((md2/(n-k))^(3/2)))
                Srho[k] <- Srho[k]/smax
            }
        }
    }
    out <- new("Srho")
    out@.Data      <- 0.5*Srho
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    out@lags       <- 1:lag.max
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}
 ## ***************************************************************************************************

Srho.biv.R <- function(x,y,lag.max,stationary=TRUE,plot=FALSE,nor=FALSE){
    ## An implementation of the bivariate entropy-based dependence measure S[rho]
    ## proposed by Granger et al. (2004) Journal of Time Series Analysis.
    ## Only deals with integer/categorical time series, for continuous processes
    ##  non-parametric density estimation is required.
    ## Probabilities are estimated through relative frequencies

    ## R VERSION OF Srho.biv.F

    ## PARAMETERS:
    ##   x,y:        integer/categorical series
    ##   lag.max:       number of lags  -- default: (round(n/4))
    ##   plot:       if TRUE plots Srho

    ## OUTPUT:
    ## S(k) = Srho(X_t,Y_{t+k)}

#**************************

Srho.biva <- function(tx,ty,jd,nor){
    dim.jd <- dim(jd)
    nx     <- dim.jd[1]
    ny     <- dim.jd[2]
    S <- 0;
    for(ix in 1:nx) {
        for(iy in 1:ny) {
            S <- S + (sqrt(jd[ix,iy])-sqrt(tx[ix]*ty[iy]))^2;
        }
    }
    if(nor){
        smax <- max(1 - sum(tx**1.5),1 - sum(ty**1.5))
        S    <- S/smax
    }
    return(S/2);
}
#**************************


    ## Check for equal length **********
    nx <- length(x);
    ny <- length(y)
    n <- min(nx,ny)
    x <- x[1:n];
    y <- y[1:n];
    ## Check for equal length **********

    Srho <- rep(0,(2*lag.max+1))

    jd  <- prop.table(table(x,y));
    tx  <- margin.table(jd,1)
    ix <- length(tx)
    ty  <- margin.table(jd,2)
    iy <- length(ty)
    Srho[lag.max+1] <- Srho.biva(tx,ty,jd,nor);
    if(stationary==TRUE) {
        for(k in 1:lag.max) {
            x1  <- x[1:(n-k)];
            y1  <- y[(k+1):n];
            dum  <- prop.table(table(x1,y1));
            jd <- matrix(0,ix,iy)
            dimnames(jd) <- list(names(tx),names(ty))
            jd[dimnames(dum)$x1,dimnames(dum)$y1] <- dum
            Srho[lag.max+1+k] <- Srho.biva(tx,ty,jd,nor);

            x1  <- x[(k+1):n];
            y1  <- y[1:(n-k)];
            jd  <- prop.table(table(x1,y1));
            Srho[lag.max+1-k] <- Srho.biva(tx,ty,jd,nor);
        }
    }
    else {
        for(k in 1:lag.max) {
            x1  <- x[1:(n-k)];
            y1  <- y[(k+1):n];
            jd  <- prop.table(table(x1,y1));
            tx  <- margin.table(jd,1)
            ty  <- margin.table(jd,2)
            Srho[lag.max+1+k] <- Srho.biva(tx,ty,jd,nor);

            x1  <- x[(k+1):n];
            y1  <- y[1:(n-k)];
            jd  <- prop.table(table(x1,y1));
            tx  <- margin.table(jd,1)
            ty  <- margin.table(jd,2)
            Srho[lag.max+1-k] <- Srho.biva(tx,ty,jd,nor);
        }
    }
    out <- new("Srho")
    out@.Data      <- Srho
    out@lags       <- -lag.max:lag.max
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

 ## ***************************************************************************************************

Srho.biv.F <- function(x,y,lag.max,stationary=TRUE,plot=FALSE,nor=FALSE){

    ## An implementation of the entropy-based dependence measure S[rho]
    ## proposed by Granger et al. (2004) Journal of Time Series Analysis.
    ## Only deals with integer/categorical time series, for continuous processes
    ## non-parametric density estimation is required.
    ## Probabilities are estimated through relative frequencies

    ## FORTRAN VERSION OF Srho.biv.R

    ## PARAMETERS:
    ##   x,y:        categorical series
    ##   lag.max:       number of lags  -- default: (round(n/4))
    ##   plot:       if TRUE plots Srho

    ## OUTPUT:
    ## S(k) = Srho(X_t,Y_{t+k)}

    ## Check for equal length **********
    nx <- length(x);
    ny <- length(y)
    n <- min(nx,ny)
    x <- x[1:n];
    y <- y[1:n];
    ## Check for equal length **********

    if(stationary==TRUE) {
        SS <- .Fortran("ssbiv2",as.integer(x),as.integer(y),as.integer(n),as.integer(lag.max),S=double((2*lag.max+1)),as.integer(nor),PACKAGE="tseriesEntropy")$S
    }
    else {
        SS <- .Fortran("ssbiv", as.integer(x),as.integer(y),as.integer(n),as.integer(lag.max),S=double((2*lag.max+1)),as.integer(nor),PACKAGE="tseriesEntropy")$S
    }
    out <- new("Srho")
    out@.Data      <-  SS
    out@lags       <- -lag.max:lag.max
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

 ## ***************************************************************************************************
