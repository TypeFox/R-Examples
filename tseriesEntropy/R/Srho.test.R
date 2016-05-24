
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

 setClass("Srho.test", contains="Srho",
         representation(call="call",call.h="call",quantiles ="matrix", test.type ="character", significant.lags = "list", p.value = "numeric"))


## ***************************************************************************************************
 setMethod("plot" , signature(x = "Srho.test",y = "missing"),
     function(x, y, type = "s", xlab = "lag", ylab = "S", ylim= c(0,max(max(x@.Data),max(x@quantiles[,1]),
     max(x@quantiles[,2]))), main = NULL, col=1, mai=c(.85,.75,.1,.1), lwd=1.5, lty.l=c(2,2), col.l=c(3,4), grid=TRUE, ...){
            par(mai= mai);
            plot(x@lags,x@.Data,type=type,col=col,xlab=xlab,ylab=ylab,ylim=ylim,lwd=lwd,...);
            abline(h=0,lty=2,col=1);
            lines(x@lags,x@quantiles[,1],lty=lty.l[1],col=col.l[1],type="s",lwd=lwd);
            lines(x@lags,x@quantiles[,2],lty=lty.l[2],col=col.l[2],type="s",lwd=lwd);
            if(grid){grid()}
      }
)
## ***************************************************************************************************
 setMethod ("show" , "Srho.test",
    function(object){
    out <- object@.Data
    names(out) <- object@lags;
    n <- length(out);
    lag.min <- object@lags[1]
    lag.max <- object@lags[n]
    cat (" -------------------------------------------------------------------------- \n")
    cat (" Srho test for", object@test.type, "on lags", lag.min, "to", lag.max, "\n")
    cat (" -------------------------------------------------------------------------- \n")
    cat (" Call: \n")
    print(object@call);
    cat (" -------------------------------------------------------------------------- \n")
    cat (" Stationary version      :" , object@stationary , "\n")
    cat (" Significant.lags:", "\n")
    print(object@significant.lags)
    cat (" -------------------------------------------------------------------------- \n")
    cat (" p-values:", "\n")
    print(object@p.value)
    cat (" -------------------------------------------------------------------------- \n")
    }
)

## ***************************************************************************************************

Srho.test <- function(x,y,lag.max,B=1000,stationary=TRUE,plot=TRUE,quant=c(0.95,0.99),nor=FALSE){

    n <- length(x);
    if(!(is.integer(x)||is.factor(x))) stop('input series must be either integer or categorical valued')
    if (missing(lag.max)) lag.max = round(n/4);
    if((lag.max >= n)||(lag.max < 2)) stop('incorrect number of lags')

    if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
    if(length(quant)==1){
        if(quant==0.99){
            quant <- c(0.95,quant)
        }else{
            quant <- c(quant,0.99)
        }
    }
    if (missing(y)){
        return(Srho.test.uni(x,lag.max=lag.max,B=B,stationary=stationary,plot=plot,quant=quant,nor=nor))
    } else {
        return(Srho.test.biv(x,y,lag.max=lag.max,B=B,stationary=stationary,plot=plot,quant=quant,nor=nor))
    }
}


## ***************************************************************************************************

Srho.test.uni <- function(x,lag.max,B=1000,stationary=TRUE,plot=TRUE,quant,nor=FALSE){

    ## Bootstrap/permutation test of serial dependence based on Srho (FORTRAN version)

    ## PARAMETERS:
    ##   x:          integer or categorical series
    ##   lag.max:    number of lags  -- default: round(n/4)
    ##   B:          number of Monte Carlo replications for the bands -- default: 1000
    ##   plot:       if TRUE plots Srho together with the confidence bands at 95% and 99%
    ## Simone Giannerini APRIL 2007

    n     <- length(x);
    M     <- matrix(0,lag.max,B);
    S.b   <- .Fortran("ssunib",as.integer(x),as.integer(n),as.integer(lag.max),B=as.integer(B), S=as.double(rep(0,lag.max)),
    M=matrix(as.double(0),lag.max,B),as.integer(stationary),as.integer(nor),PACKAGE="tseriesEntropy")

    S.x   <- S.b$S;
    M     <- S.b$M;
    M.95  <- apply(M,1,quantile,probs=c(quant[1]));
    M.99  <- apply(M,1,quantile,probs=c(quant[2]));
    names(S.x) <- 1:lag.max
    ind95 <- which(S.x>=M.95);
    ind99 <- which(S.x>=M.99);
    out   <- new("Srho.test")
    out@.Data      <- S.x
    out@lags       <- 1:lag.max
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    out@test.type  <- "serial dependence"
    out@quantiles  <- cbind(M.95,M.99)
    q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
    colnames(out@quantiles)     <- q.names
    rownames(out@quantiles)     <- 1:lag.max
    out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
    names(out@significant.lags) <- q.names
    out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
    names(out@p.value)          <- 1:lag.max
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

## ***************************************************************************************************

Srho.test.biv <- function(x,y,lag.max,B=1000,stationary=TRUE,plot=TRUE,quant,nor=FALSE){

   ## Bivariate bootstrap/permutation test of serial dependence based on the Cross entropy Srho (FORTRAN version)

    ## PARAMETERS:
    ##   x:          integer or categorical series
    ##   y:          integer or categorical series
    ##   lag.max:       number of lags  -- default: (round(n/4)
    ##   B:          number of Monte Carlo replications for the bands -- default: 2000
    ##   plot:       if TRUE plots Srho together with the confidence bands at 95% and 99%

    ## OUTPUT:
    ## S(k) = Srho(X_t,Y_{t+k)}

    ## Simone Giannerini December 2007

    ## Check for equal length **********
    nx <- length(x);
    ny <- length(y)
    if(nx!=ny) warning('The two series have different lengths')
    n  <- min(nx,ny)
    x  <- x[1:n];
    y  <- y[1:n];
    ## Check for equal length **********

    M   <- matrix(0,(2*lag.max+1),B);storage.mode(M)<-"double"
    S.b <- .Fortran("ssbivb",as.integer(x),as.integer(y),as.integer(n),as.integer(lag.max),B=as.integer(B),
        S=double((2*lag.max+1)),M=M,as.integer(stationary),as.integer(nor),PACKAGE="tseriesEntropy");
    #SUBROUTINE SSBIVB(X,Y,N,nlag,B,S,M,STAT,nor)
    S.x   <- S.b$S;
    M     <- S.b$M;
    M.95  <- apply(M,1,quantile,probs=c(quant[1]));
    M.99  <- apply(M,1,quantile,probs=c(quant[2]));
    names(S.x) <- -lag.max:lag.max
    ind95 <- which(S.x>=M.95);
    ind99 <- which(S.x>=M.99);
    out   <- new("Srho.test")
    out@.Data <- S.x
    out@lags  <- -lag.max:lag.max
    out@stationary <- stationary
    out@data.type  <- "integer-categorical"
    out@test.type  <- "cross dependence"
    out@quantiles  <- cbind(M.95,M.99)
    q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
    colnames(out@quantiles)     <- q.names
    rownames(out@quantiles)     <- -lag.max:lag.max
    out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
    names(out@significant.lags) <- q.names
    out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
    names(out@p.value)          <- -lag.max:lag.max
    if(nor){out@notes <- "normalized"}
    if (plot) {
        plot(out)
        return(invisible(out))
    }
    else return(out)
}

 ## ***************************************************************************************************
