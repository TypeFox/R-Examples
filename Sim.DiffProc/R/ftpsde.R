## Tue Jan 12 15:47:04 2016
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
##### fptsde1d

fptsde1d <- function(N, ...)  UseMethod("fptsde1d")

fptsde1d.default <- function(N =1000,M=100,x0=0,t0=0,T=1,Dt,boundary,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x' ")
    if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
    if (length(all.vars(boundary)) > 1 )  stop("boundary must depend on 't' or constant")
    if (length(all.vars(boundary)) == 1 && all.vars(boundary) != "t" ) stop("boundary must depend on 't' or constant")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0 ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt  <- (T - t0)/N 
    res <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type,method,...)
    R   <- data.frame(res$X)
    Bn  <- function(t)  eval(boundary)+0*t 
    F   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R[,i])) )
      if (x0 > Bn(t0)){
           v1  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(res))))==0, NA,min(which(R[,i] <= Bn(time(res))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v1[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v1[i]],tol= .Machine$double.eps)$root))
      }else if (x0 < Bn(t0)){
           v2  <- sapply(1:length(F),function(i) ifelse(length(which(R[,i] <= Bn(time(res))))==(N+1), NA,min(which(R[,i] >= Bn(time(res))))))
           fpt <- sapply(1:length(F),function(i) ifelse(is.na(v2[i]),NA,uniroot(f = function(x) (F[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v2[i]],tol= .Machine$double.eps)$root))
      }else{
           fpt <- rep(0,M)
                     }
    structure(list(SDE=res,boundary=boundary[[1]],fpt=fpt),class="fptsde1d")
}

###

bconfint.fptsde1d <- function(x,level=0.95,...)
                 {
    class(x) <-"fptsde1d"
    x <- x$fpt
    bconfint(x,level=level,...)
}

skewness.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$fpt
    skewness(x) 
}

kurtosis.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$fpt
    kurtosis(x) 
}

median.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$fpt
    median(x,na.rm = TRUE) 
}

mean.fptsde1d <- function(x,...)
                    {
    class(x) <- "fptsde1d"
    x <- x$fpt
    mean(x,na.rm = TRUE) 
}

quantile.fptsde1d <- function(x,...)
             {
    class(x) <- "fptsde1d"
    x <- x$fpt
    quantile(x,na.rm = TRUE,...)
}

moment.fptsde1d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde1d"
    x <- x$fpt
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) moment(x,order=order[j]))))
    rownames(Mom) <- paste("fpt",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.fptsde1d <- function(object, ...)
           {
    x <- object	   
    class(x) <- "fptsde1d"
    S <- function(t) eval(x$boundary)
    if (x$SDE$x0 > S(x$SDE$t0)){
    cat("\n\t Monte-Carlo Statistics for the F.P.T of X(t)","\n", "| T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) <= ", deparse(x$boundary),"}","\n",
        sep="")
    }else{
    cat("\n\t Monte-Carlo Statistics for the F.P.T of X(t)","\n", "| T(S,X) = inf{t >= ",x$SDE$t0 ," : X(t) >= ", deparse(x$boundary),"}","\n",
        sep="")
    }
    res <- data.frame(matrix(c(length(which(is.na(x$fpt))),sprintf("%f",mean(x$fpt,na.rm = TRUE)),sprintf("%f",var(x$fpt,na.rm = TRUE)),
                               sprintf("%f",median(x$fpt,na.rm = TRUE)),sprintf("%f",quantile(x$fpt,0.25,na.rm = TRUE)),
							   sprintf("%f",quantile(x$fpt,0.75,na.rm = TRUE)),sprintf("%f",skewness(x$fpt)),sprintf("%f",kurtosis(x$fpt)),
							   sprintf("%f",moment(x$fpt,order=2)),sprintf("%f",moment(x$fpt,order=3)),sprintf("%f",moment(x$fpt,order=4)),
							   sprintf("%f",moment(x$fpt,order=5)),sprintf("%f",bconfint(x$fpt)[1]),sprintf("%f",bconfint(x$fpt)[2])),
                               ncol=1))
    # rownames(res) <- paste(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c(""),sep="")
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("fpt(x)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(x)
}

plot.fptsde1d <- function(x,...) .plot.fptsde1d(x,...)

################################################################################
################################################################################
#####
##### fptsde2d

fptsde2d <- function(N, ...)  UseMethod("fptsde2d")


fptsde2d.default <- function(N =1000,M=100,x0=0,y0=0,t0=0,T=1,Dt,boundary,driftx,diffx,
                     drifty,diffy,alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler",
                     "milstein","predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
    if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
    if (length(all.vars(boundary)) > 1 )  stop("boundary must depend on 't' or constant")
    if (length(all.vars(boundary)) == 1 && all.vars(boundary) != "t" ) stop("boundary must depend on 't' or constant")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N 
    res <- snssde2d(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)
    R1   <- data.frame(res$X)
    R2   <- data.frame(res$Y)
    Bn  <- function(t)  eval(boundary)+0*t 
    F1   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R1[,i])) )
    F2   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R2[,i])) )
    if (x0 > Bn(t0)){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(res))))==0, NA,min(which(R1[,i] <= Bn(time(res))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (x0 < Bn(t0)){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(res))))==(N+1), NA,min(which(R1[,i] >= Bn(time(res))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,M)
                     }
    if (y0 > Bn(t0)){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(res))))==0, NA,min(which(R2[,i] <= Bn(time(res))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (y0 < Bn(t0)){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(res))))==(N+1), NA,min(which(R2[,i] >= Bn(time(res))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,M)
                     }
    structure(list(SDE=res,boundary=boundary[[1]],fptx=fptx,fpty=fpty),class="fptsde2d")
}

###

bconfint.fptsde2d <- function(x,level=0.95,...)
             {
    class(x) <- "fptsde2d"
    Bcon <- t(data.frame(do.call("cbind",lapply(3:4,function(i) bconfint(x[[i]][!is.na(x[[i]])],level=level)))))
    rownames(Bcon) <- paste(c("fpt(x)","fpt(y)"),sep="")
    return(Bcon)
}

skewness.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Skew <- data.frame(do.call("cbind",lapply(3:4,function(i) skewness(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Skew) <-  paste(c("fpt(x)","fpt(y)"),sep="")
    return(Skew)
}

kurtosis.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Kurt <- data.frame(do.call("cbind",lapply(3:4,function(i) kurtosis(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Kurt) <- paste(c("fpt(x)","fpt(y)"),sep="")
    return(Kurt)
}

median.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Med <- data.frame(do.call("cbind",lapply(3:4,function(i) median(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Med) <- paste(c("fpt(x)","fpt(y)"),sep="")
    return(Med)
}

mean.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Mean <- data.frame(do.call("cbind",lapply(3:4,function(i) mean(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Mean) <- paste(c("fpt(x)","fpt(y)"),sep="")
    return(Mean)
}

quantile.fptsde2d <- function(x,...)
             {
    class(x) <- "fptsde2d"
    Qun <- t(data.frame(do.call("cbind",lapply(3:4,function(i) quantile(x[[i]][!is.na(x[[i]])],...)))))
    rownames(Qun) <- paste(c("fpt(x)","fpt(y)"),sep="")
    return(Qun)
}

moment.fptsde2d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde2d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(3:4,function(i) moment(x[[i]][!is.na(x[[i]])],order=order[j])))))
    rownames(Mom) <- paste(c("fpt(x)","fpt(y)"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}


summary.fptsde2d <- function(object, ...)
           {  
    class(object) <- "fptsde2d"
    S <- function(t) eval(object$boundary)
    if (object$SDE$x0 > S(object$SDE$t0)){
    fpt_x <- paste("T(S,X) = inf{t >= ",object$SDE$t0 ," : X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S,X) = inf{t >= ",object$SDE$t0 ," : X(t) >= ", deparse(object$boundary),"}")
    }
    if (object$SDE$y0 > S(object$SDE$t0)){
    fpt_y <- paste("T(S,Y) = inf{t >= ",object$SDE$t0 ," : Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S,Y) = inf{t >= ",object$SDE$t0 ," : Y(t) >= ", deparse(object$boundary),"}")
    }
    cat("\n\t Monte-Carlo Statistics for the F.P.T of (X(t),Y(t))","\n",
        "| ",fpt_x,"\n",
        "| ",fpt_y,"\n\n",
       sep="")
    x <- object$fptx##[!is.na(X$fptx)]
    y <- object$fpty##[!is.na(X$fpty)]
    res <- data.frame(matrix(c(length(which(is.na(x))),sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2]),
                               length(which(is.na(y))),sprintf("%f",mean(y,na.rm = TRUE)),sprintf("%f",var(y,na.rm = TRUE)),sprintf("%f",median(y,na.rm = TRUE)),
                               sprintf("%f",quantile(y,0.25,na.rm = TRUE)),sprintf("%f",quantile(y,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(y)),sprintf("%f",kurtosis(y)),sprintf("%f",moment(y,order=2)),sprintf("%f",moment(y,order=3)),
                               sprintf("%f",moment(y,order=4)),sprintf("%f",moment(y,order=5)),sprintf("%f",bconfint(y)[1]),sprintf("%f",bconfint(y)[2])),
                               ncol=2))
    # rownames(res) <- paste(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("fpt(x)","fpt(y)"),sep="")
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("fpt(x)","fpt(y)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

plot.fptsde2d <- function(x,...) .plot.fptsde2d(x,...)
	

################################################################################
################################################################################
#####
##### fptsde3d

fptsde3d <- function(N, ...)  UseMethod("fptsde3d")

fptsde3d.default <- function(N =1000,M=100,x0=0,y0=0,z0=0,t0=0,T=1,Dt,boundary,driftx,diffx,drifty,diffy,
                     driftz,diffz,alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler",
                     "milstein","predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (any(!is.numeric(x0) || !is.numeric(y0) || !is.numeric(z0))) stop("'x0','y0' and 'z0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(boundary))) stop(" must be expressions of a constant or time-dependent boundary ")
    if (length(all.vars(boundary)) > 1 )  stop("boundary must depend on 't' or constant")
    if (length(all.vars(boundary)) == 1 && all.vars(boundary) != "t" ) stop("boundary must depend on 't' or constant")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if(any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if(any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (any(t0 < 0 || T < 0 || T <= t0) ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N 
    res <- snssde3d(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)
    R1   <- data.frame(res$X)
    R2   <- data.frame(res$Y)
    R3   <- data.frame(res$Z)	
    Bn  <- function(t)  eval(boundary)+0*t 
    F1   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R1[,i])) )
    F2   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R2[,i])) )
    F3   <- lapply(1:M,function(i) approxfun(time(res),as.vector(R3[,i])) )
    if (x0 > Bn(t0)){
           v11  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(res))))==0, NA,min(which(R1[,i] <= Bn(time(res))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v11[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v11[i]],tol= .Machine$double.eps)$root))
      }else if (x0 < Bn(t0)){
           v12  <- sapply(1:length(F1),function(i) ifelse(length(which(R1[,i] <= Bn(time(res))))==(N+1), NA,min(which(R1[,i] >= Bn(time(res))))))
           fptx <- sapply(1:length(F1),function(i) ifelse(is.na(v12[i]),NA,uniroot(f = function(x) (F1[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v12[i]],tol= .Machine$double.eps)$root))
      }else{
           fptx <- rep(0,M)
                     }
    if (y0 > Bn(t0)){
           v21  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(res))))==0, NA,min(which(R2[,i] <= Bn(time(res))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v21[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v21[i]],tol= .Machine$double.eps)$root))
      }else if (y0 < Bn(t0)){
           v22  <- sapply(1:length(F2),function(i) ifelse(length(which(R2[,i] <= Bn(time(res))))==(N+1), NA,min(which(R2[,i] >= Bn(time(res))))))
           fpty <- sapply(1:length(F2),function(i) ifelse(is.na(v22[i]),NA,uniroot(f = function(x) (F2[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v22[i]],tol= .Machine$double.eps)$root))
      }else{
           fpty <- rep(0,M)
                     }
    if (z0 > Bn(t0)){
           v31  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(res))))==0, NA,min(which(R3[,i] <= Bn(time(res))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v31[i]),NA,uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v31[i]],tol= .Machine$double.eps)$root))
      }else if (z0 < Bn(t0)){
           v32  <- sapply(1:length(F3),function(i) ifelse(length(which(R3[,i] <= Bn(time(res))))==(N+1), NA,min(which(R3[,i] >= Bn(time(res))))))
           fptz <- sapply(1:length(F3),function(i) ifelse(is.na(v32[i]),NA,uniroot(f = function(x) (F3[[i]](x)- Bn(x)),lower=t0,upper=time(res)[v32[i]],tol= .Machine$double.eps)$root))
      }else{
           fptz <- rep(0,M)
                     }
    structure(list(SDE=res,boundary=boundary[[1]],fptx=fptx,fpty=fpty,fptz=fptz),class="fptsde3d")
}

###

bconfint.fptsde3d <- function(x,level=0.95,...)
             {
    class(x) <- "fptsde3d"
    Bcon <- t(data.frame(do.call("cbind",lapply(3:5,function(i) bconfint(x[[i]][!is.na(x[[i]])],level=level)))))
    rownames(Bcon) <- paste(c("fptx","fpty","fptz"),sep="")
    return(Bcon)
}

skewness.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Skew <- data.frame(do.call("cbind",lapply(3:5,function(i) skewness(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Skew) <- paste(c("fptx","fpty","fptz"),sep="")
    return(Skew)
}

kurtosis.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Kurt <- data.frame(do.call("cbind",lapply(3:5,function(i) kurtosis(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Kurt) <-paste(c("fptx","fpty","fptz"),sep="")
    return(Kurt)
}

median.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Med <- data.frame(do.call("cbind",lapply(3:5,function(i) median(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Med) <-paste(c("fptx","fpty","fptz"),sep="")
    return(Med)
}

mean.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Mean <- data.frame(do.call("cbind",lapply(3:5,function(i) mean(x[[i]][!is.na(x[[i]])]))),rownames = "")
    names(Mean) <-paste(c("fptx","fpty","fptz"),sep="")
    return(Mean)
}

quantile.fptsde3d <- function(x,...)
             {
    class(x) <- "fptsde3d"
    Qun <- t(data.frame(do.call("cbind",lapply(3:5,function(i) quantile(x[[i]][!is.na(x[[i]])],...)))))
    rownames(Qun) <-paste(c("fptx","fpty","fptz"),sep="")
    return(Qun)
}

moment.fptsde3d <- function(x,order = 2,...)
             {
    class(x) <- "fptsde3d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(3:5,function(i) moment(x[[i]][!is.na(x[[i]])],order=order[j])))))
    rownames(Mom) <-paste(c("fptx","fpty","fptz"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.fptsde3d <- function(object, ...)
           {  
    class(object) <- "fptsde3d"
    S <- function(t) eval(object$boundary)
    if (object$SDE$x0 > S(object$SDE$t0)){
    fpt_x <- paste("T(S,X) = inf{t >= ",object$SDE$t0 ," : X(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_x <- paste("T(S,X) = inf{t >= ",object$SDE$t0 ," : X(t) >= ", deparse(object$boundary),"}")
    }
    if (object$SDE$y0 > S(object$SDE$t0)){
    fpt_y <- paste("T(S,Y) = inf{t >= ",object$SDE$t0 ," : Y(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_y <- paste("T(S,Y) = inf{t >= ",object$SDE$t0 ," : Y(t) >= ", deparse(object$boundary),"}")
    }
    if (object$SDE$z0 > S(object$SDE$t0)){
    fpt_z <- paste("T(S,Z) = inf{t >= ",object$SDE$t0 ," : Z(t) <= ", deparse(object$boundary),"}")
    }else{
    fpt_z <- paste("T(S,Z) = inf{t >= ",object$SDE$t0 ," : Z(t) >= ", deparse(object$boundary),"}")
    }
    cat("\n\tMonte-Carlo Statistics for the F.P.T of (X(t),Y(t),Z(t))","\n",
        "| ",fpt_x,"\n",
        "| ",fpt_y,"\n",
        "| ",fpt_z,"\n\n",
       sep="")
    x <- object$fptx##[!is.na(X$fptx)]
    y <- object$fpty##[!is.na(X$fpty)]
    z <- object$fptz##[!is.na(X$fpty)]
    res <- data.frame(matrix(c(length(which(is.na(x))),sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2]),
                               length(which(is.na(y))),sprintf("%f",mean(y,na.rm = TRUE)),sprintf("%f",var(y,na.rm = TRUE)),sprintf("%f",median(y,na.rm = TRUE)),
                               sprintf("%f",quantile(y,0.25,na.rm = TRUE)),sprintf("%f",quantile(y,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(y)),sprintf("%f",kurtosis(y)),sprintf("%f",moment(y,order=2)),sprintf("%f",moment(y,order=3)),
                               sprintf("%f",moment(y,order=4)),sprintf("%f",moment(y,order=5)),sprintf("%f",bconfint(y)[1]),sprintf("%f",bconfint(y)[2]),
							   length(which(is.na(z))),sprintf("%f",mean(z,na.rm = TRUE)),sprintf("%f",var(z,na.rm = TRUE)),sprintf("%f",median(z,na.rm = TRUE)),
                               sprintf("%f",quantile(z,0.25,na.rm = TRUE)),sprintf("%f",quantile(z,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(z)),sprintf("%f",kurtosis(z)),sprintf("%f",moment(z,order=2)),sprintf("%f",moment(z,order=3)),
                               sprintf("%f",moment(z,order=4)),sprintf("%f",moment(z,order=5)),sprintf("%f",bconfint(z)[1]),sprintf("%f",bconfint(z)[2])),
                               ncol=3))
    # rownames(res) <- paste(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("fpt(x)","fpt(y)","fpt(z)"),sep="")
	dimnames(res) <- list(c("NA's","Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("fpt(x)","fpt(y)","fpt(z)"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

plot.fptsde3d <- function(x,...) .plot.fptsde3d(x,...)

