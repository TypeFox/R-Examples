## Tue Sep 09 13:54:04 2014
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
##### bridgesde1d

bridgesde1d <- function(N, ...)  UseMethod("bridgesde1d")

bridgesde1d.default <- function(N =1000,M=1,x0=0,y=0,t0=0,T=1,Dt,drift,diffusion,
                              alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || !is.numeric(y)) stop("'x0' and 'y' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N
	X1 <- snssde1d(N,M,x0,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X
    if (M > 1){X2 <- apply(data.frame(snssde1d(N,M,x0=y,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X),2,rev)}else{X2 <- rev(snssde1d(N,M,x0=y,t0,T,Dt,drift,diffusion,alpha,mu,type, method,...)$X)}        
    G <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        G[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        G[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        G <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        G <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   G[which(G==0)]=NA
   name <- "X"
   name <- if(M > 1) paste("X",1:length(which(!is.na(G))),sep="")
   if (M == 1){ if (is.na(G) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{X <- ts(c(X1[1:G],X2[-(1:G)]),start=t0,deltat=deltat(X1),names=name)}
   }else if (M > 1){ if(length(which(is.na(G))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(G))) == 1){X <- ts(c(X1[,-which(is.na(G))][1:G[which(!is.na(G))]],X2[,-which(is.na(G))][-(1:G[which(!is.na(G))])]),start=t0,deltat=deltat(X1),names=name)
                 }else if (length(which( is.na(G))) == 0){X <- ts(sapply(1:length(which(!is.na(G))), function(j) c(X1[,j][1:G[!is.na(G)][j]],X2[,j][-(1:G[!is.na(G)][j])])),start=t0,deltat=deltat(X1),names=name)
                 }else{ X1 <- X1[,-c(which(is.na(G)))]
				        X2 <- X2[,-c(which(is.na(G)))]
                        X <- ts(sapply(1:length(which(!is.na(G))), function(j) c(X1[,j][1:G[!is.na(G)][j]],X2[,j][-(1:G[!is.na(G)][j])])),
                                 start=t0,deltat=deltat(X1),names=name)
                       }
   }
    structure(list(X=X,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T,C=G),class="bridgesde1d")
}

###

print.bridgesde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde1d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 2"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
	Dr <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$drift), ignore.case = F,fixed = T)
    DD <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$diffusion), ignore.case = F,fixed = T)
    if(x$type=="ito"){
    cat("Ito Bridges Sde 1D:","\n",        
        "\t| dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| C  = ",format(length(which(!is.na(x$C))),digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Final value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Bridges Sde 1D:","\n",
        "\t| dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| C  = ",format(length(which(!is.na(x$C))),digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Final value","\t\t| y = ",format(x$y,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

time.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
    as.vector(time(x$X))
}

# summary.bridgesde1d  <- function(object, ...)
           # {   
    # class(object) <- "bridgesde1d"
    # cat("\n\tMonte-Carlo Statistics for X(t) at final time T = ",object$T,"\n\n",
        # sep="")
    # if (length(which(!is.na(object$C))) == 1 ){
    # x <- object$X[which(time(object)==object$T)]}else{
    # x <- object$X[which(time(object)==object$T),]}
    # res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               # sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               # sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               # sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2])),
                               # ncol=1))
    # row.names(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          # "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  # "Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("X"),sep="")
    # print(res, quote = FALSE, right = TRUE,...)
    # invisible(object)
# }


mean.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    rowMeans(Y,na.rm = TRUE,...)
}

skewness.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Skew <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    row.names(Skew) <- paste("X(t=",time(x),")",sep="")
    names(Skew) <- paste(c("Skewness"),sep="")
    return(Skew)
}

kurtosis.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    kurt <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    row.names(kurt) <- paste("X(t=",time(x),")",sep="")
    names(kurt) <- paste(c("Kurtosis"),sep="")
    return(kurt)
}

median.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Med <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    row.names(Med) <- paste("X(t=",time(x),")",sep="")
    names(Med) <- paste(c("Median"),sep="")
    return(Med)
}

quantile.bridgesde1d <- function(x,...)
                    {
    class(x) <- "bridgesde1d"
    if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Qun <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    row.names(Qun) <- paste("X(t=",time(x),")",sep="")
    return(Qun)
}

moment.bridgesde1d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Mom <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    row.names(Mom) <- paste("X(t=",time(x),")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

bconfint.bridgesde1d <- function(x,level = 0.95,...)
                    {
    class(x) <- "bridgesde1d"
	if (length(which(!is.na(x$C))) == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    conf <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(conf) <- paste("X(t=",time(x),")",sep="")
    names(conf) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(conf)
}

##
## Plot

plot.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
    plot(X,...)
}

lines.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
	if (length(which(!is.na(x$C))) >=2){
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}}else{
	lines(time(x),X,...)}
}


points.bridgesde1d <- function(x,...)
                 {
    class(x) <- "bridgesde1d"
    X <- x$X
	if (length(which(!is.na(x$C))) >=2){
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}}else{
	points(time(x),X,...)}
}


#####
##### bridgesde2d

bridgesde2d <- function(N, ...)  UseMethod("bridgesde2d")

bridgesde2d.default <- function(N =1000,M=1,x0=c(0,0),y=c(1,1),t0=0,T=1,Dt,driftx,diffx,drifty,
                              diffy,alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 2) stop("'x0' must be numeric, and length(x0) = 2 ")
    if (!is.numeric(y)  || length(y) != 2) stop("'y' must be numeric, and length(y) = 2 ")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x' and 'y'")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if (t0 < 0 || T < 0) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N
    X1 <- snssde2d(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$X
    Y1 <- snssde2d(N,M,x0=x0[1],y0=x0[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$Y
    if (M > 1){X2 <- apply(data.frame(snssde2d(N,M,x0=y[1],y0=y[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$X),2,rev)
	           Y2 <- apply(data.frame(snssde2d(N,M,x0=y[1],y0=y[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$Y),2,rev)
    }else{X2 <- rev(snssde2d(N,M,x0=y[1],y0=y[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$X)
	      Y2 <- rev(snssde2d(N,M,x0=y[1],y0=y[2],t0,T,Dt,driftx,diffx,drifty,diffy,alpha,mu,type, method,...)$Y)
	    }        
    Gx = Gy <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        Gx[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        Gx[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        Gx <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        Gx <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Y1[1,M] >= Y2[1,M]){
                    if (!all(Y1[,j] > Y2[,j]))
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) <= 0)) - 1
             }else{ if (!all(Y1[,j] < Y2[,j])) 
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Y1[1] >= Y2[1]){
                    if (!all(Y1 > Y2))
                        Gy <- min(which((Y1-Y2) <= 0)) - 1
            }else{  if (!all(Y1 < Y2)) 
                        Gy <- min(which((Y1-Y2) >= 0)) - 1
						}       
   }
   Gx[which(Gx==0)]=NA
   Gy[which(Gy==0)]=NA
   namex <- "X"
   namey <- "Y"
   namex <- if(M > 1) paste("X",1:length(which(!is.na(Gx))),sep="")
   namey <- if(M > 1) paste("Y",1:length(which(!is.na(Gy))),sep="")
   if (M == 1){ if (is.na(Gx) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{X <- ts(c(X1[1:Gx],X2[-(1:Gx)]),start=t0,deltat=deltat(X1),names=namex)}
   }else if (M > 1){ if(length(which(is.na(Gx))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gx))) == 1){X <- ts(c(X1[,-which(is.na(Gx))][1:Gx[which(!is.na(Gx))]],X2[,-which(is.na(Gx))][-(1:Gx[which(!is.na(Gx))])]),start=t0,deltat=deltat(X1),names=namex)
                 }else if (length(which( is.na(Gx))) == 0){X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),start=t0,deltat=deltat(X1),names=namex)
                 }else{ X1 <- X1[,-c(which(is.na(Gx)))]
				        X2 <- X2[,-c(which(is.na(Gx)))]
                        X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),
                                 start=t0,deltat=deltat(X1),names=namex)
                       }
   }
   if (M == 1){ if (is.na(Gy) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Y <- ts(c(Y1[1:Gy],Y2[-(1:Gy)]),start=t0,deltat=deltat(Y1),names=namey)}
   }else if (M > 1){ if(length(which(is.na(Gy))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gy))) == 1){Y <- ts(c(Y1[,-which(is.na(Gy))][1:Gy[which(!is.na(Gy))]],Y2[,-which(is.na(Gy))][-(1:Gy[which(!is.na(Gy))])]),start=t0,deltat=deltat(Y1),names=namey)
                 }else if (length(which( is.na(Gy))) == 0){Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),start=t0,deltat=deltat(Y1),names=namey)
                 }else{ Y1 <- Y1[,-c(which(is.na(Gy)))]
				        Y2 <- Y2[,-c(which(is.na(Gy)))]
                        Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),
                                 start=t0,deltat=deltat(Y1),names=namey)
                       }
   }
    structure(list(X=X,Y=Y, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],type=type,method=method, 
                   x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T,Cx=Gx,Cy=Gy),class="bridgesde2d")
}

###

print.bridgesde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde2d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 1.5"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 2"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
	Drx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	Dry <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDy <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = as.expression(x$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    if(x$type=="ito"){
    cat("Ito Bridges Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| (Cx,Cy) = c","(",format(length(which(!is.na(x$Cx))),digits=digits),",",format(length(which(!is.na(x$Cy))),digits=digits),")",".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridges Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| (Cx,Cy) = c","(",format(length(which(!is.na(x$Cx))),digits=digits),",",format(length(which(!is.na(x$Cy))),digits=digits),")",".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.bridgesde2d <- function(x,...) .plot.bridgesde2d(x,...)

lines.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}				  
    for (i in 1:length(which(!is.na(x$Cx)))){lines(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){lines(time(x),Y[,i],...)}
}

points.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}				  
    for (i in 1:length(which(!is.na(x$Cx)))){points(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){points(time(x),Y[,i],...)}
}

plot2d.bridgesde2d <- function(x,...) .plot2d.bridgesde2d(x,...)

lines2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    X <- X[,1]
    Y <- Y[,1]
    lines2d(X,Y,...)
}

points2d.bridgesde2d <- function(x,...)
                 {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    X <- X[,1]
    Y <- Y[,1]
    points2d(X,Y,...)
}

##
## summary

# summary.bridgesde2d <- function(object,...)
                    # {
	# x <- object				
    # class(x) <- "bridgesde2d"
    # summary(x$XY,...)
# }

mean.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    return(list(X=rowMeans(X,na.rm = TRUE,...),Y=rowMeans(Y,na.rm = TRUE,...)))
}

skewness.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    Skewx <- data.frame(sapply(1:(x$N+1),function(i) skewness(X[i,]) ))
    row.names(Skewx) <- paste("X(t=",time(x),")",sep="")
    names(Skewx) <- paste(c("Skewness"),sep="")
    Skewy <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    row.names(Skewy) <- paste("Y(t=",time(x),")",sep="")
    names(Skewy) <- paste(c("Skewness"),sep="")
    return(list(X=Skewx,Y=Skewy))
}

kurtosis.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    kurtx <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(X[i,]) ))
    row.names(kurtx) <- paste("X(t=",time(x),")",sep="")
    names(kurtx) <- paste(c("Kurtosis"),sep="")
    kurty <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    row.names(kurty) <- paste("Y(t=",time(x),")",sep="")
    names(kurty) <- paste(c("Kurtosis"),sep="")
    return(list(X=kurtx,Y=kurty))
}

median.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    Medx <- data.frame(sapply(1:(x$N+1),function(i) median(X[i,],na.rm = TRUE) ))
    row.names(Medx) <- paste("X(t=",time(x),")",sep="")
    names(Medx) <- paste(c("Median"),sep="")
    Medy <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    row.names(Medy) <- paste("Y(t=",time(x),")",sep="")
    names(Medy) <- paste(c("Median"),sep="")
    return(list(X=Medx,Y=Medy))
}

quantile.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    Qunx <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(X[i,],na.rm = TRUE,...) ))))
    row.names(Qunx) <- paste("X(t=",time(x),")",sep="")
    Quny <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    row.names(Quny) <- paste("Y(t=",time(x),")",sep="")
    return(list(X=Qunx,Y=Quny))
}

moment.bridgesde2d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    Momx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(X[i,],order=order[j],...)))))
    row.names(Momx) <- paste("X(t=",time(x),")",sep="")
    names(Momx) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    row.names(Momy) <- paste("Y(t=",time(x),")",sep="")
    names(Momy) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(list(X=Momx,Y=Momy))
}

bconfint.bridgesde2d <- function(x,level = 0.95,...)
                    {
    class(x) <- "bridgesde2d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    confx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(confx) <- paste("X(t=",time(x),")",sep="")
    names(confx) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    confy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(confy) <- paste("Y(t=",time(x),")",sep="")
    names(confy) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(list(X=confx,Y=confy))
}


time.bridgesde2d <- function(x,...)
                    {
    class(x) <- "bridgesde2d"
    as.vector(time(x$X))
}

#####
##### bridgesde3d

bridgesde3d <- function(N, ...)  UseMethod("bridgesde3d")

bridgesde3d.default <- function(N =1000,M=1,x0=c(0,0,0),y=c(1,-1,2),t0=0,T=1,Dt,driftx,diffx,drifty,
                              diffy,driftz,diffz,alpha=0.5,mu=0.5,type=c("ito","str"), method=c(
                              "euler","milstein","predcorr","smilstein","taylor",
                              "heun","rk1","rk2","rk3"),...)
                     {
    if (!is.numeric(x0) || length(x0) != 3) stop("'x0' must be numeric, and length(x0) = 3 ")
    if (!is.numeric(y)  || length(y) != 3) stop("'y' must be numeric, and length(y) = 3 ")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (missing(type)) type <- "ito"
    method <- match.arg(method)
    if (method =="predcorr"){
    if (any(alpha > 1 || alpha < 0)) stop("please use '0 <= alpha <= 1' ")
    if (any(mu > 1 || mu < 0))       stop("please use '0 <= mu <= 1' ")
                            }
    if ( t0 < 0 || T < 0 ) stop(" please use positive times! (0 <= t0 < T) ")
    if (missing(Dt)) {
        t <- seq(t0, T, length = N + 1)
    } else {
        t <- c(t0, t0 + cumsum(rep(Dt, N)))
        T <- t[N + 1]
    }
    Dt <- (T - t0)/N	
	res1 <- snssde3d(N,M=1,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)
	res2 <- 
	
	X1 <- snssde3d(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$X
    Y1 <- snssde3d(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Y
	Z1 <- snssde3d(N,M,x0=x0[1],y0=x0[2],z0=x0[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Z
    if (M > 1){X2 <- apply(data.frame(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$X),2,rev)
	           Y2 <- apply(data.frame(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Y),2,rev)
			   Z2 <- apply(data.frame(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Z),2,rev)
    }else{X2 <- rev(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$X)
	      Y2 <- rev(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Y)
	      Z2 <- rev(snssde3d(N,M,x0=y[1],y0=y[2],z0=y[3],t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,alpha,mu,type, method,...)$Z)		  
	    }        
    Gx = Gy = Gz <- rep(NA,M)
    if (M > 1){
        for(j in 1:M){
              if (X1[1,M] >= X2[1,M]){
                    if (!all(X1[,j] > X2[,j]))
                        Gx[j] <- min(which((X1[,j]-X2[,j]) <= 0)) - 1
             }else{ if (!all(X1[,j] < X2[,j])) 
                        Gx[j] <- min(which((X1[,j]-X2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (X1[1] >= X2[1]){
                    if (!all(X1 > X2))
                        Gx <- min(which((X1-X2) <= 0)) - 1
            }else{  if (!all(X1 < X2)) 
                        Gx <- min(which((X1-X2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Y1[1,M] >= Y2[1,M]){
                    if (!all(Y1[,j] > Y2[,j]))
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) <= 0)) - 1
             }else{ if (!all(Y1[,j] < Y2[,j])) 
                        Gy[j] <- min(which((Y1[,j]-Y2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Y1[1] >= Y2[1]){
                    if (!all(Y1 > Y2))
                        Gy <- min(which((Y1-Y2) <= 0)) - 1
            }else{  if (!all(Y1 < Y2)) 
                        Gy <- min(which((Y1-Y2) >= 0)) - 1
						}       
   }
   if (M > 1){
        for(j in 1:M){
              if (Z1[1,M] >= Z2[1,M]){
                    if (!all(Z1[,j] > Z2[,j]))
                        Gz[j] <- min(which((Z1[,j]-Z2[,j]) <= 0)) - 1
             }else{ if (!all(Z1[,j] < Z2[,j])) 
                        Gz[j] <- min(which((Z1[,j]-Z2[,j]) >= 0)) - 1
						}
                   }
    }else{
             if (Z1[1] >= Z2[1]){
                    if (!all(Z1 > Z2))
                        Gz <- min(which((Z1-Z2) <= 0)) - 1
            }else{  if (!all(Z1 < Z2)) 
                        Gz <- min(which((Z1-Z2) >= 0)) - 1
						}       
   }   
   Gx[which(Gx==0)]=NA
   Gy[which(Gy==0)]=NA
   Gz[which(Gz==0)]=NA
   namex <- "X"
   namey <- "Y"
   namey <- "Z"
   namex <- if(M > 1) paste("X",1:length(which(!is.na(Gx))),sep="")
   namey <- if(M > 1) paste("Y",1:length(which(!is.na(Gy))),sep="")
   namez <- if(M > 1) paste("Z",1:length(which(!is.na(Gz))),sep="")
   if (M == 1){ if (is.na(Gx) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{X <- ts(c(X1[1:Gx],X2[-(1:Gx)]),start=t0,deltat=deltat(X1),names=namex)}
   }else if (M > 1){ if(length(which(is.na(Gx))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gx))) == 1){X <- ts(c(X1[,-which(is.na(Gx))][1:Gx[which(!is.na(Gx))]],X2[,-which(is.na(Gx))][-(1:Gx[which(!is.na(Gx))])]),start=t0,deltat=deltat(X1),names=namex)
                 }else if (length(which( is.na(Gx))) == 0){X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),start=t0,deltat=deltat(X1),names=namex)
                 }else{ X1 <- X1[,-c(which(is.na(Gx)))]
				        X2 <- X2[,-c(which(is.na(Gx)))]
                        X <- ts(sapply(1:length(which(!is.na(Gx))), function(j) c(X1[,j][1:Gx[!is.na(Gx)][j]],X2[,j][-(1:Gx[!is.na(Gx)][j])])),
                                 start=t0,deltat=deltat(X1),names=namex)
                       }
   }
   if (M == 1){ if (is.na(Gy) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Y <- ts(c(Y1[1:Gy],Y2[-(1:Gy)]),start=t0,deltat=deltat(Y1),names=namey)}
   }else if (M > 1){ if(length(which(is.na(Gy))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gy))) == 1){Y <- ts(c(Y1[,-which(is.na(Gy))][1:Gy[which(!is.na(Gy))]],Y2[,-which(is.na(Gy))][-(1:Gy[which(!is.na(Gy))])]),start=t0,deltat=deltat(Y1),names=namey)
                 }else if (length(which( is.na(Gy))) == 0){Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),start=t0,deltat=deltat(Y1),names=namey)
                 }else{ Y1 <- Y1[,-c(which(is.na(Gy)))]
				        Y2 <- Y2[,-c(which(is.na(Gy)))]
                        Y <- ts(sapply(1:length(which(!is.na(Gy))), function(j) c(Y1[,j][1:Gy[!is.na(Gy)][j]],Y2[,j][-(1:Gy[!is.na(Gy)][j])])),
                                 start=t0,deltat=deltat(Y1),names=namey)
                       }
   }
   if (M == 1){ if (is.na(Gz) ){stop( "A crossing has been no realized,trying again (Repeat)..." )}else{Z <- ts(c(Z1[1:Gz],Z2[-(1:Gz)]),start=t0,deltat=deltat(Z1),names=namez)}
   }else if (M > 1){ if(length(which(is.na(Gz))) == M){stop( "A crossing has been no realized,trying again (Repeat)..." )
                 }else if (length(which(!is.na(Gz))) == 1){Z <- ts(c(Z1[,-which(is.na(Gz))][1:Gz[which(!is.na(Gz))]],Z2[,-which(is.na(Gz))][-(1:Gz[which(!is.na(Gz))])]),start=t0,deltat=deltat(Z1),names=namez)
                 }else if (length(which( is.na(Gz))) == 0){Z <- ts(sapply(1:length(which(!is.na(Gz))), function(j) c(Z1[,j][1:Gz[!is.na(Gz)][j]],Z2[,j][-(1:Gz[!is.na(Gz)][j])])),start=t0,deltat=deltat(Z1),names=namez)
                 }else{ Z1 <- Z1[,-c(which(is.na(Gz)))]
				        Z2 <- Z2[,-c(which(is.na(Gz)))]
                        Z <- ts(sapply(1:length(which(!is.na(Gz))), function(j) c(Z1[,j][1:Gz[!is.na(Gz)][j]],Z2[,j][-(1:Gz[!is.na(Gz)][j])])),
                                 start=t0,deltat=deltat(Z1),names=namez)
                       }
   }
   structure(list(X=X,Y=Y,Z=Z, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],driftz=driftz[[1]],diffz=diffz[[1]],
                   type=type,method=method, x0=x0,y=y, N=N,Dt=Dt,t0=t0,T=T,Cx=Gx,Cy=Gy,Cz=Gz),class="bridgesde3d")
}

###

print.bridgesde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "bridgesde3d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 2"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
    Drx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDx <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    Dry <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDy <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	Drz <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$driftz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDz <- gsub(pattern = 'x', replacement = 'X(t)', x = gsub(pattern = 'y', replacement = 'Y(t)', x = gsub(pattern = 'z', replacement = 'Z(t)', x = as.expression(x$diffz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    if(x$type=="ito"){
    cat("Ito Bridges Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| (Cx,Cy,Cz) = c","(",format(length(which(!is.na(x$Cx))),digits=digits),",",format(length(which(!is.na(x$Cy))),digits=digits),",",format(length(which(!is.na(x$Cz))),digits=digits),")",".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Bridges Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Crossing realized","\t| (Cx,Cy,Cz) = c","(",format(length(which(!is.na(x$Cx))),digits=digits),",",format(length(which(!is.na(x$Cy))),digits=digits),",",format(length(which(!is.na(x$Cz))),digits=digits),")",".","\n",
        "\t| Initial values","\t| x0 = c","(",format(x$x0[1],digits=digits),",",format(x$x0[2],digits=digits),",",format(x$x0[3],digits=digits),")",".","\n",
        "\t| Final values","\t\t| y  = c","(",format(x$y[1],digits=digits),",",format(x$y[2],digits=digits),",",format(x$y[3],digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## summary

# summary.bridgesde3d <- function(object,...)
                    # {
	# x <- object				
    # class(x) <- "bridgesde3d"
    # summary(x$XYZ,...)
# }

time.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    as.vector(time(x$X))
}

mean.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    return(list(X=rowMeans(X,na.rm = TRUE,...),Y=rowMeans(Y,na.rm = TRUE,...),Z=rowMeans(Z,na.rm = TRUE,...)))
}

skewness.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    Skewx <- data.frame(sapply(1:(x$N+1),function(i) skewness(X[i,]) ))
    row.names(Skewx) <- paste("X(t=",time(x),")",sep="")
    names(Skewx) <- paste(c("Skewness"),sep="")
    Skewy <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    row.names(Skewy) <- paste("Y(t=",time(x),")",sep="")
    names(Skewy) <- paste(c("Skewness"),sep="")
    Skewz <- data.frame(sapply(1:(x$N+1),function(i) skewness(Z[i,]) ))
    row.names(Skewz) <- paste("Z(t=",time(x),")",sep="")
    names(Skewz) <- paste(c("Skewness"),sep="")	
    return(list(X=Skewx,Y=Skewy,Z=Skewz))
}

kurtosis.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    kurtx <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(X[i,]) ))
    row.names(kurtx) <- paste("X(t=",time(x),")",sep="")
    names(kurtx) <- paste(c("Kurtosis"),sep="")
    kurty <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    row.names(kurty) <- paste("Y(t=",time(x),")",sep="")
    names(kurty) <- paste(c("Kurtosis"),sep="")
    kurtz <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Z[i,]) ))
    row.names(kurtz) <- paste("Z(t=",time(x),")",sep="")
    names(kurtz) <- paste(c("Kurtosis"),sep="")
    return(list(X=kurtx,Y=kurty,Z=kurtz))
}

median.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    Medx <- data.frame(sapply(1:(x$N+1),function(i) median(X[i,],na.rm = TRUE) ))
    row.names(Medx) <- paste("X(t=",time(x),")",sep="")
    names(Medx) <- paste(c("Median"),sep="")
    Medy <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    row.names(Medy) <- paste("Y(t=",time(x),")",sep="")
    names(Medy) <- paste(c("Median"),sep="")
    Medz <- data.frame(sapply(1:(x$N+1),function(i) median(Z[i,],na.rm = TRUE) ))
    row.names(Medz) <- paste("Z(t=",time(x),")",sep="")
    names(Medz) <- paste(c("Median"),sep="")
    return(list(X=Medx,Y=Medy,Z=Medz))
}

quantile.bridgesde3d <- function(x,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    Qunx <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(X[i,],na.rm = TRUE,...) ))))
    row.names(Qunx) <- paste("X(t=",time(x),")",sep="")
    Quny <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    row.names(Quny) <- paste("Y(t=",time(x),")",sep="")
    Qunz <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Z[i,],na.rm = TRUE,...) ))))
    row.names(Qunz) <- paste("Z(t=",time(x),")",sep="")
    return(list(X=Qunx,Y=Quny,Z=Qunz))
}

moment.bridgesde3d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    Momx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(X[i,],order=order[j],...)))))
    row.names(Momx) <- paste("X(t=",time(x),")",sep="")
    names(Momx) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    row.names(Momy) <- paste("Y(t=",time(x),")",sep="")
    names(Momy) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momz <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Z[i,],order=order[j],...)))))
    row.names(Momz) <- paste("Z(t=",time(x),")",sep="")
    names(Momz) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(list(X=Momx,Y=Momy,Z=Momz))
}

bconfint.bridgesde3d <- function(x,level = 0.95,...)
                    {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}
    confx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(confx) <- paste("X(t=",time(x),")",sep="")
    names(confx) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    confy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(confy) <- paste("Y(t=",time(x),")",sep="")
    names(confy) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
	confz <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Z[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    row.names(confz) <- paste("Z(t=",time(x),")",sep="")
    names(confz) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(list(X=confx,Y=confy,Z=confz))
}

##
## Plot

plot.bridgesde3d <- function(x,...) .plot.bridgesde3d(x,...)

lines.bridgesde3d <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}		
    for (i in 1:length(which(!is.na(x$Cx)))){lines(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){lines(time(x),Y[,i],...)} 
    for (i in 1:length(which(!is.na(x$Cz)))){lines(time(x),Z[,i],...)} 	
}

points.bridgesde3d <- function(x,...)
                 {
    class(x) <- "bridgesde3d"
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}		
    for (i in 1:length(which(!is.na(x$Cx)))){points(time(x),X[,i],...)}
    for (i in 1:length(which(!is.na(x$Cy)))){points(time(x),Y[,i],...)} 
    for (i in 1:length(which(!is.na(x$Cz)))){points(time(x),Z[,i],...)} 
}


plot3D.bridgesde3d <- function(x,display = c("persp","rgl"),...)
                 {
	class(x) <- "bridgesde3d"		 
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    if (length(which(!is.na(x$Cz))) == 1){Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{Z = x$Z}	
    plot3D(X[,1],Y[,1],Z[,1],display,...)
}

