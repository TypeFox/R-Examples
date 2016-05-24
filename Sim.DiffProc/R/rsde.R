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
##### rsde1d

rsde1d <- function(N, ...)  UseMethod("rsde1d")

rsde1d.default <- function(N =1000,M=100,x0=0,t0=0,T=1,Dt,tau=0.5,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")
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
    if (any(T < tau || t0 > tau) )  stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 
    res <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)
    if (M == 1){  X = matrix(res$X,nrow=length(res$X),ncol=1)}else{X = res$X}
    X   <- as.vector(X[which(time(res)==tau),])
    if (length(X) == 0){
	if (M==1){ F   <- lapply(1:M,function(i) approxfun(time(res),res$X))}else{
               F   <- lapply(1:M,function(i) approxfun(time(res),res$X[,i]))}
               X   <- sapply(1:length(F),function(i) F[[i]](tau)) 
    }
    structure(list(SDE=res,tau=tau,x=X),class="rsde1d")
}

###

bconfint.rsde1d <- function(x,level=0.95,...)
                 {
     class(x) <-"rsde1d"
     bconfint(x$x,level=level,...)
}

skewness.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    skewness(x$x) 
}

kurtosis.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    kurtosis(x$x) 
}

median.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    median(x$x,na.rm = TRUE) 
}

mean.rsde1d <- function(x,...)
                    {
    class(x) <- "rsde1d"
    mean(x$x,na.rm = TRUE) 
}

quantile.rsde1d <- function(x,...)
             {
    class(x) <- "rsde1d"
    quantile(x$x,na.rm = TRUE,...)
}

moment.rsde1d <- function(x,order = 2,...)
             {
    class(x) <- "rsde1d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) moment(x$x,order=order[j]))))
    rownames(Mom) <- paste("X(t=",x$t,")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde1d <- function(object, ...)
           {	   
    class(object) <- "rsde1d"
    cat("\n\tMonte-Carlo Statistics for X(t) at t = ",object$tau,"\n",
        sep="")
    x <- object$x		
    res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),
                               sprintf("%f",median(x,na.rm = TRUE)),sprintf("%f",quantile(x,0.25,na.rm = TRUE)),
							   sprintf("%f",quantile(x,0.75,na.rm = TRUE)),sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),
							   sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),sprintf("%f",moment(x,order=4)),
							   sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2])),
                               ncol=1))
    # rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("x"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("x"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

plot.rsde1d <- function(x,...) .plot.rsde1d(x,...)

#####
##### rsde2d

rsde2d <- function(N, ...)  UseMethod("rsde2d")

rsde2d.default <- function(N =1000,M=100,x0=0,y0=0,t0=0,T=1,Dt,tau=0.5,driftx,diffx,drifty,diffy,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
                    { 
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T)))  stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
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
    if (any(T < tau || t0 > tau) )   stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 		
	res <- snssde2d(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)		
    if (M == 1){  X = matrix(res$X,nrow=length(res$X),ncol=1)
	              Y = matrix(res$Y,nrow=length(res$Y),ncol=1)}else{
				  X = res$X
				  Y = res$Y}
    X   <- as.vector(X[which(time(res)==tau),])
    Y   <- as.vector(Y[which(time(res)==tau),])
    if (length(X) == 0){
	if (M==1){ Fx   <- lapply(1:M,function(i) approxfun(time(res),res$X))}else{
               Fx   <- lapply(1:M,function(i) approxfun(time(res),res$X[,i]))}
               X   <- sapply(1:length(Fx),function(i) Fx[[i]](tau)) 
    }
    if (length(Y) == 0){
	if (M==1){ Fy   <- lapply(1:M,function(i) approxfun(time(res),res$Y))}else{
               Fy   <- lapply(1:M,function(i) approxfun(time(res),res$Y[,i]))}
               Y   <- sapply(1:length(Fy),function(i) Fy[[i]](tau)) 
    }
    structure(list(SDE=res,tau=tau,x=X,y=Y),class="rsde2d")
}

##
##

bconfint.rsde2d <- function(x,level=0.95,...)
             {
    class(x) <- "rsde2d"
    Bcon <- t(data.frame(do.call("cbind",lapply(3:4,function(i) bconfint(x[i][[1]],level=level)))))
    rownames(Bcon) <- paste(c("x","y"),sep="")
    return(Bcon)
}

skewness.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Skew <- data.frame(do.call("cbind",lapply(3:4,function(i) skewness(x[i][[1]]))),rownames = "")
    names(Skew) <- paste(c("x","y"),sep="")
    return(Skew)
}

kurtosis.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Kurt <- data.frame(do.call("cbind",lapply(3:4,function(i) kurtosis(x[i][[1]]))),rownames = "")
    names(Kurt) <- paste(c("x","y"),sep="")
    return(Kurt)
}

median.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Med <- data.frame(do.call("cbind",lapply(3:4,function(i) median(x[i][[1]]))),rownames = "")
    names(Med) <- paste(c("x","y"),sep="")
    return(Med)
}

mean.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Mean <- data.frame(do.call("cbind",lapply(3:4,function(i) mean(x[i][[1]]))),rownames = "")
    names(Mean) <- paste(c("x","y"),sep="")
    return(Mean)
}

quantile.rsde2d <- function(x,...)
             {
    class(x) <- "rsde2d"
    Qun <- t(data.frame(do.call("cbind",lapply(3:4,function(i) quantile(x[i][[1]],...)))))
    rownames(Qun) <- paste(c("x","y"),sep="")
    return(Qun)
}

moment.rsde2d <- function(x,order = 2,...)
             {
    class(x) <- "rsde2d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(3:4,function(i) moment(x[i][[1]],order=order[j])))))
    rownames(Mom) <- paste(c("x","y"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde2d <- function(object, ...)
           {   
    class(object) <- "rsde2d"
    cat("\n    Monte-Carlo Statistics for X(t) and Y(t) at t = ",object$tau,"\n\n",
        sep="")
    x <- object$x##[!is.na(X$fptx)]
    y <- object$y##[!is.na(X$fpty)]
    res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2]),
                               sprintf("%f",mean(y,na.rm = TRUE)),sprintf("%f",var(y,na.rm = TRUE)),sprintf("%f",median(y,na.rm = TRUE)),
                               sprintf("%f",quantile(y,0.25,na.rm = TRUE)),sprintf("%f",quantile(y,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(y)),sprintf("%f",kurtosis(y)),sprintf("%f",moment(y,order=2)),sprintf("%f",moment(y,order=3)),
                               sprintf("%f",moment(y,order=4)),sprintf("%f",moment(y,order=5)),sprintf("%f",bconfint(y)[1]),sprintf("%f",bconfint(y)[2])),
                               ncol=2))
    # rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("x","y"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("x","y"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

plot.rsde2d <- function(x,...) .plot.rsde2d(x,...)

#####
##### rsde3d

rsde3d <- function(N, ...)  UseMethod("rsde3d")

rsde3d.default <- function(N =1000,M=100,x0=0,y0=0,z0=0,t0=0,T=1,Dt,tau=0.5,driftx,diffx,drifty,diffy,driftz,diffz,
                             alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler","milstein",
                             "predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
                    { 
    if (any(!is.numeric(x0) || !is.numeric(y0) || !is.numeric(z0))) stop("'x0','y0' and 'z0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
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
    if (any(T < tau || t0 > tau) )   stop( " please use 't0 <= tau <= T'")
    Dt <- (T - t0)/N 
	res <- snssde3d(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)		
    if (M == 1){  X = matrix(res$X,nrow=length(res$X),ncol=1)
	              Y = matrix(res$Y,nrow=length(res$Y),ncol=1)
				  Z = matrix(res$Z,nrow=length(res$Z),ncol=1)}else{
				  X = res$X
				  Y = res$Y
				  Z = res$Z}
    X   <- as.vector(X[which(time(res)==tau),])
    Y   <- as.vector(Y[which(time(res)==tau),])
    Z   <- as.vector(Z[which(time(res)==tau),])
    if (length(X) == 0){
	if (M==1){ Fx   <- lapply(1:M,function(i) approxfun(time(res),res$X))}else{
               Fx   <- lapply(1:M,function(i) approxfun(time(res),res$X[,i]))}
               X   <- sapply(1:length(Fx),function(i) Fx[[i]](tau)) 
    }
    if (length(Y) == 0){
	if (M==1){ Fy   <- lapply(1:M,function(i) approxfun(time(res),res$Y))}else{
               Fy   <- lapply(1:M,function(i) approxfun(time(res),res$Y[,i]))}
               Y   <- sapply(1:length(Fy),function(i) Fy[[i]](tau)) 
    }
    if (length(Z) == 0){
	if (M==1){ Fz   <- lapply(1:M,function(i) approxfun(time(res),res$Z))}else{
               Fz   <- lapply(1:M,function(i) approxfun(time(res),res$Z[,i]))}
               Z   <- sapply(1:length(Fz),function(i) Fz[[i]](tau)) 
    }
    structure(list(SDE=res,tau=tau,x=X,y=Y,z=Z),class="rsde3d")
}

###

bconfint.rsde3d <- function(x,level=0.95,...)
             {
    class(x) <- "rsde3d"
    Bcon <- t(data.frame(do.call("cbind",lapply(3:5,function(i) bconfint(x[i][[1]],level=level)))))
    rownames(Bcon) <- paste(c("x","y","z"),sep="")
    return(Bcon)
}

skewness.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Skew <- data.frame(do.call("cbind",lapply(3:5,function(i) skewness(x[i][[1]]))),rownames = "")
    names(Skew) <- paste(c("x","y","z"),sep="")
    return(Skew)
}

kurtosis.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Kurt <- data.frame(do.call("cbind",lapply(3:5,function(i) kurtosis(x[i][[1]]))),rownames = "")
    names(Kurt) <- paste(c("x","y","z"),sep="")
    return(Kurt)
}

median.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Med <- data.frame(do.call("cbind",lapply(3:5,function(i) median(x[i][[1]]))),rownames = "")
    names(Med) <- paste(c("x","y","z"),sep="")
    return(Med)
}

mean.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Mean <- data.frame(do.call("cbind",lapply(3:5,function(i) mean(x[i][[1]]))),rownames = "")
    names(Mean) <- paste(c("x","y","z"),sep="")
    return(Mean)
}

quantile.rsde3d <- function(x,...)
             {
    class(x) <- "rsde3d"
    Qun <- t(data.frame(do.call("cbind",lapply(3:5,function(i) quantile(x[i][[1]],...)))))
    rownames(Qun) <- paste(c("x","y","z"),sep="")
    return(Qun)
}

moment.rsde3d <- function(x,order = 2,...)
             {
    class(x) <- "rsde3d"
    Mom <- data.frame(do.call("cbind",lapply(1:length(order), function(j) sapply(3:5,function(i) moment(x[i][[1]],order=order[j])))))
    rownames(Mom) <- paste(c("x","y","z"),sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

summary.rsde3d <- function(object, ...)
           {	   
    class(object) <- "rsde3d"
    cat("\n Monte-Carlo Statistics for X(t), Y(t) and Z(t) at t = ",object$tau,"\n\n",
        sep="")
    x <- object$x##[!is.na(X$fptx)]
    y <- object$y##[!is.na(X$fpty)]
    z <- object$z##[!is.na(X$fpty)]
    res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2]),
                               sprintf("%f",mean(y,na.rm = TRUE)),sprintf("%f",var(y,na.rm = TRUE)),sprintf("%f",median(y,na.rm = TRUE)),
                               sprintf("%f",quantile(y,0.25,na.rm = TRUE)),sprintf("%f",quantile(y,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(y)),sprintf("%f",kurtosis(y)),sprintf("%f",moment(y,order=2)),sprintf("%f",moment(y,order=3)),
                               sprintf("%f",moment(y,order=4)),sprintf("%f",moment(y,order=5)),sprintf("%f",bconfint(y)[1]),sprintf("%f",bconfint(y)[2]),
							   sprintf("%f",mean(z,na.rm = TRUE)),sprintf("%f",var(z,na.rm = TRUE)),sprintf("%f",median(z,na.rm = TRUE)),
                               sprintf("%f",quantile(z,0.25,na.rm = TRUE)),sprintf("%f",quantile(z,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(z)),sprintf("%f",kurtosis(z)),sprintf("%f",moment(z,order=2)),sprintf("%f",moment(z,order=3)),
                               sprintf("%f",moment(z,order=4)),sprintf("%f",moment(z,order=5)),sprintf("%f",bconfint(z)[1]),sprintf("%f",bconfint(z)[2])),
                               ncol=3))
    # rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile",
                              # "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              # "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("x","y","z"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile",
                              "Skewness","Kurtosis","Moment of order 2","Moment of order 3",
                              "Moment of order 4","Moment of order 5","Bound conf Inf (95%)","Bound conf Sup (95%)"),c("x","y","z"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

plot.rsde3d <- function(x,...) .plot.rsde3d(x,...)

