## Tue Jan 12 15:37:25 2016
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
##### snssde1D

snssde1d <- function(N, ...)  UseMethod("snssde1d")

snssde1d.default <- function(N =1000,M=1,x0=0,t0=0,T=1,Dt,drift,diffusion,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (!is.numeric(x0)) stop("'x0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (missing(drift))     drift     <- expression(0)
    if (missing(diffusion)) diffusion <- expression(1)
    if (!is.expression(drift) || !is.expression(diffusion)) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't' and 'x'")
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
    if (method=="euler")         {res <- .Euler1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="predcorr") {res <- .PredCorr1D(N,M,x0,t0,T,Dt,alpha,mu,drift,diffusion,type)}
    else if (method=="milstein") {res <- .Milstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="smilstein"){res <- .SMilstein1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="taylor")   {res <- .STS1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="heun")     {res <- .Heun1D(N,M,x0,t0,T,Dt,drift,diffusion,type)}
    else if (method=="rk1")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=1)}
    else if (method=="rk2")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=2)}
    else if (method=="rk3")      {res <- .RK1D(N,M,x0,t0,T,Dt,drift,diffusion,type,order=3)}
    structure(list(X=res$X,drift=drift[[1]], diffusion=diffusion[[1]],type=type,method=method, 
                   x0=x0, N=N, M=M,Dt=Dt,t0=t0,T=T),class="snssde1d")
}

###

print.snssde1d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde1d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme"}
    else if (x$method=="taylor")   {sch <- "Ito-Taylor scheme of order 1.5"}
    else if (x$method=="heun")     {sch <- "Heun scheme of order 2"}
    else if (x$method=="rk1")      {sch <- "Runge-Kutta method of order 1"}
    else if (x$method=="rk2")      {sch <- "Runge-Kutta method of order 2"}
    else if (x$method=="rk3")      {sch <- "Runge-Kutta method of order 3"}
	Dr <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$drift), ignore.case = F,fixed = T)
    DD <- gsub(pattern = 'x', replacement = 'X(t)', x = as.expression(x$diffusion), ignore.case = F,fixed = T)
    if(x$type=="ito"){
    cat("Ito Sde 1D:","\n",        
        "\t| dX(t) = ", Dr," * dt + ", DD," * dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",       
        sep="")}else{
    cat("Stratonovich Sde 1D:","\n",
        "\t| dX(t) = ", Dr," * dt + ", DD," o dW(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial value","\t\t| x0 = ",format(x$x0,digits=digits),".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    plot(X,...)
}

lines.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
	if (x$M >=2){
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}}else{
	lines(time(x),X,...)}
}


points.snssde1d <- function(x,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
	if (x$M >=2){
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}}else{
	points(time(x),X,...)}
}

add.mean <- function(x,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 {
    class(x) <- "snssde1d"
    X <- x$X
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(col)) {col = 2}
    if (is.null(cex)) {cex = 0.8}
	if (x$M >=2){
    lines(time(x),rowMeans(X,na.rm = TRUE),lwd=lwd,lty=lty,col=col,...)}else{
	lines(time(x),X,lwd=lwd,lty=lty,col=col,...)}
    legend("topright",c("mean path"),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
}


add.bconfint.snssde1d <- function(x,level=0.95,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 {
    class(x) <- "snssde1d"
    if (is.null(lty)) {lty = 1}
    if (is.null(lwd)) {lwd = 1}
    if (is.null(col)) {col = 4}
    if (is.null(cex)) {cex = 0.8}
    lines(time(x),bconfint(x,level)[,1],lwd=lwd,lty=lty,col=col,...)
    lines(time(x),bconfint(x,level)[,2],lwd=lwd,lty=lty,col=col,...)
    legend("topleft",c(paste("bound of",level*100,"% confidence")),inset = .01,lty=lty,col=col,lwd=lwd,cex=cex,...)
}

##
## summary

summary.snssde1d  <- function(object, ...)
           {   
    class(object) <- "snssde1d"
    cat("\n\tMonte-Carlo Statistics for X(t) at final time T = ",object$T,"\n\n",
        sep="")
    if (object$M == 1 ){
    x <- object$X[which(time(object)==object$T)]}else{
    x <- object$X[which(time(object)==object$T),]}
    res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2])),
                               ncol=1))
    #rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	#                          "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
	#						  "Bound conf Sup (95%)"),sep="")
    #names(res) <- paste(c("X"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  "Bound conf Sup (95%)"),c("X"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}


mean.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    rowMeans(Y,na.rm = TRUE,...)
}

skewness.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Skew <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    rownames(Skew) <- paste("X(t=",time(x),")",sep="")
    names(Skew) <- paste(c("Skewness"),sep="")
    return(Skew)
}

kurtosis.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    kurt <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    rownames(kurt) <- paste("X(t=",time(x),")",sep="")
    names(kurt) <- paste(c("Kurtosis"),sep="")
    return(kurt)
}

median.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Med <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    rownames(Med) <- paste("X(t=",time(x),")",sep="")
    names(Med) <- paste(c("Median"),sep="")
    return(Med)
}

quantile.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Qun <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    rownames(Qun) <- paste("X(t=",time(x),")",sep="")
    return(Qun)
}

moment.snssde1d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    Mom <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    rownames(Mom) <- paste("X(t=",time(x),")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

bconfint.snssde1d <- function(x,level = 0.95,...)
                    {
    class(x) <- "snssde1d"
	if (x$M == 1){Y = matrix(x$X,nrow=length(x$X),ncol=1)}else{Y = x$X}
    conf <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(conf) <- paste("X(t=",time(x),")",sep="")
    names(conf) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(conf)
}

#var.x <- function(x,...)
#                    {
#    class(x) <- "snssde1d" 
#    apply(x$X,1,var)
#}

time.snssde1d <- function(x,...)
                    {
    class(x) <- "snssde1d"
    as.vector(time(x$X))
}


################################################################################
################################################################################
##### snssde2d

snssde2d <- function(N, ...)  UseMethod("snssde2d")

snssde2d.default <- function(N =1000,M=1,x0=0,y0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,alpha=0.5,mu=0.5,
                     type=c("ito","str"), method=c("euler","milstein","predcorr",
                     "smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || !is.numeric(y0))) stop("'x0' and 'y0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
	if (missing(driftx) & missing(drifty) ) driftx = drifty <- expression(0)
    if (missing(diffx) & missing(diffy))    diffx = diffy   <- expression(1)
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
    if (method=="euler")         {res <- .Euler2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="predcorr") {res <- .PredCorr2D(N,M,x0,y0,t0,T,Dt,alpha,mu,driftx,diffx,drifty,diffy,type)}
    else if (method=="milstein") {res <- .Milstein2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="smilstein"){res <- .SMilstein2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="taylor")   {res <- .STS2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="heun")     {res <- .Heun2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type)}
    else if (method=="rk1")      {res <- .RK2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=1)}
    else if (method=="rk2")      {res <- .RK2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=2)}
    else if (method=="rk3")      {res <- .RK2D(N,M,x0,y0,t0,T,Dt,driftx,diffx,drifty,diffy,type,order=3)}
    structure(list(X=res$X,Y=res$Y, driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],type=type,method=method, 
                   x0=x0,y0=y0, N=N,M=M,Dt=Dt,t0=t0,T=T),class="snssde2d")
}

###

print.snssde2d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde2d"
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
    cat("Ito Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 2D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde2d <- function(x,...) .plot.snssde2d(x,...)

lines.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}    
    for (i in 1:x$M){
    lines(time(x),X[,i],...)
    lines(time(x),Y[,i],...)}
}

points.snssde2d <- function(x,...)
                 {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y} 
    for (i in 1:x$M){
    points(time(x),X[,i],...)
    points(time(x),Y[,i],...)}
}

plot2d.snssde2d <- function(x,...) .plot2d.snssde2d(x,...)

lines2d.snssde2d <- function(x,...)
        {
    class(x) <- "snssde2d"
    if (x$M == 1){
    lines(as.vector(x$X),as.vector(x$Y),...)}else{
    lines(as.vector(x$X[,1]),as.vector(x$Y[,1]),...)}
}

points2d.snssde2d <- function(x,...)
        {
    class(x) <- "snssde2d"
    if (x$M == 1){
    points(as.vector(x$X),as.vector(x$Y),...)}else{
    points(as.vector(x$X[,1]),as.vector(x$Y[,1]),...)}
}

##
## summary

summary.snssde2d  <- function(object, ...)
           {   
    class(object) <- "snssde2d"
    cat("\n\tMonte-Carlo Statistics for (X(t),Y(t)) at final time T = ",object$T,"\n\n",
        sep="")
    if (object$M == 1 ){
    x <- object$X[which(time(object)==object$T)]
    y <- object$Y[which(time(object)==object$T)]}else{
    x <- object$X[which(time(object)==object$T),]
    y <- object$Y[which(time(object)==object$T),]}
    res <- data.frame(matrix(c(sprintf("%f",mean(x,na.rm = TRUE)),sprintf("%f",var(x,na.rm = TRUE)),sprintf("%f",median(x,na.rm = TRUE)),
                               sprintf("%f",quantile(x,0.25,na.rm = TRUE)),sprintf("%f",quantile(x,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(x)),sprintf("%f",kurtosis(x)),sprintf("%f",moment(x,order=2)),sprintf("%f",moment(x,order=3)),
                               sprintf("%f",moment(x,order=4)),sprintf("%f",moment(x,order=5)),sprintf("%f",bconfint(x)[1]),sprintf("%f",bconfint(x)[2]),
                               sprintf("%f",mean(y,na.rm = TRUE)),sprintf("%f",var(y,na.rm = TRUE)),sprintf("%f",median(y,na.rm = TRUE)),
                               sprintf("%f",quantile(y,0.25,na.rm = TRUE)),sprintf("%f",quantile(y,0.75,na.rm = TRUE)),
                               sprintf("%f",skewness(y)),sprintf("%f",kurtosis(y)),sprintf("%f",moment(y,order=2)),sprintf("%f",moment(y,order=3)),
                               sprintf("%f",moment(y,order=4)),sprintf("%f",moment(y,order=5)),sprintf("%f",bconfint(y)[1]),sprintf("%f",bconfint(y)[2])),
                               ncol=2))
    # rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          # "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  # "Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("X","Y"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  "Bound conf Sup (95%)"),c("X","Y"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}

mean.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    return(list(X=rowMeans(X,na.rm = TRUE,...),Y=rowMeans(Y,na.rm = TRUE,...)))
}

skewness.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    Skewx <- data.frame(sapply(1:(x$N+1),function(i) skewness(X[i,]) ))
    rownames(Skewx) <- paste("X(t=",time(x),")",sep="")
    names(Skewx) <- paste(c("Skewness"),sep="")
    Skewy <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    rownames(Skewy) <- paste("Y(t=",time(x),")",sep="")
    names(Skewy) <- paste(c("Skewness"),sep="")
    return(list(X=Skewx,Y=Skewy))
}

kurtosis.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    kurtx <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(X[i,]) ))
    rownames(kurtx) <- paste("X(t=",time(x),")",sep="")
    names(kurtx) <- paste(c("Kurtosis"),sep="")
    kurty <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    rownames(kurty) <- paste("Y(t=",time(x),")",sep="")
    names(kurty) <- paste(c("Kurtosis"),sep="")
    return(list(X=kurtx,Y=kurty))
}

median.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    Medx <- data.frame(sapply(1:(x$N+1),function(i) median(X[i,],na.rm = TRUE) ))
    rownames(Medx) <- paste("X(t=",time(x),")",sep="")
    names(Medx) <- paste(c("Median"),sep="")
    Medy <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    rownames(Medy) <- paste("Y(t=",time(x),")",sep="")
    names(Medy) <- paste(c("Median"),sep="")
    return(list(X=Medx,Y=Medy))
}

quantile.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    Qunx <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(X[i,],na.rm = TRUE,...) ))))
    rownames(Qunx) <- paste("X(t=",time(x),")",sep="")
    Quny <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    rownames(Quny) <- paste("Y(t=",time(x),")",sep="")
    return(list(X=Qunx,Y=Quny))
}

moment.snssde2d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    Momx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(X[i,],order=order[j],...)))))
    rownames(Momx) <- paste("X(t=",time(x),")",sep="")
    names(Momx) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    rownames(Momy) <- paste("Y(t=",time(x),")",sep="")
    names(Momy) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(list(X=Momx,Y=Momy))
}

bconfint.snssde2d <- function(x,level = 0.95,...)
                    {
    class(x) <- "snssde2d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y}
    confx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(confx) <- paste("X(t=",time(x),")",sep="")
    names(confx) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    confy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(confy) <- paste("Y(t=",time(x),")",sep="")
    names(confy) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(list(X=confx,Y=confy))
}

time.snssde2d <- function(x,...)
                    {
    class(x) <- "snssde2d"
    as.vector(time(x$X))
}


################################################################################
################################################################################
##### snssde3d


snssde3d <- function(N, ...)  UseMethod("snssde3d")

snssde3d.default <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,driftz,diffz,
                             alpha=0.5,mu=0.5,type=c("ito","str"), method=c("euler","milstein",
                             "predcorr","smilstein","taylor","heun","rk1","rk2","rk3"),...)
        {
    if (any(!is.numeric(x0) || !is.numeric(y0) || !is.numeric(z0))) stop("'x0', 'y0' and 'z0' must be numeric")
    if (any(!is.numeric(t0) || !is.numeric(T))) stop(" 't0' and 'T' must be numeric")
    if (any(!is.numeric(N)  || (N - floor(N) > 0) || N <= 1)) stop(" 'N' must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0))  stop(" 'M' must be a positive integer ")
    if (missing(driftx) & missing(drifty) & missing(driftz)) driftx = drifty = driftz <- expression(0)
    if (missing(diffx) & missing(diffy) & missing(diffz))    diffx = diffy = diffz  <- expression(1)
    if (any(!is.expression(driftx) || !is.expression(diffx) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(drifty) || !is.expression(diffy) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
    if (any(!is.expression(driftz) || !is.expression(diffz) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions in 't', 'x', 'y' and 'z'")
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
    if (method=="euler")         {res <- .Euler3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="predcorr") {res <- .PredCorr3D(N,M,x0,y0,z0,t0,T,Dt,alpha,mu,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="milstein") {res <- .Milstein3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="smilstein"){res <- .SMilstein3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="taylor")   {res <- .STS3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="heun")     {res <- .Heun3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type)}
    else if (method=="rk1")      {res <- .RK3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=1)}
    else if (method=="rk2")      {res <- .RK3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=2)}
    else if (method=="rk3")      {res <- .RK3D(N,M,x0,y0,z0,t0,T,Dt,driftx,diffx,drifty,diffy,driftz,diffz,type,order=3)}
    structure(list(X=res$X,Y=res$Y,Z=res$Z,driftx=driftx[[1]], diffx=diffx[[1]],drifty=drifty[[1]], diffy=diffy[[1]],driftz=driftz[[1]], 
                   diffz=diffz[[1]],type=type,method=method,x0=x0,y0=y0,z0=z0,N=N,M=M,Dt=Dt,t0=t0,T=T),class="snssde3d")
}


###

print.snssde3d <- function(x, digits=NULL, ...)
           {
    class(x) <- "snssde3d"
    if (x$method=="euler")         {sch <- "Euler scheme of order 0.5"}
    else if (x$method=="milstein") {sch <- "Milstein scheme of order 1"}
    else if (x$method=="predcorr") {sch <- "Predictor-corrector method of order 1"}
    else if (x$method=="smilstein"){sch <- "Second Milstein scheme of order 1.5"}
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
    cat("Ito Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," * dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," * dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," * dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich Sde 3D:","\n",
        "\t| dX(t) = ", Drx," * dt + ", DDx," o dW1(t)","\n", 
        "\t| dY(t) = ", Dry," * dt + ", DDy," o dW2(t)","\n",
        "\t| dZ(t) = ", Drz," * dt + ", DDz," o dW3(t)","\n",
        "Method:","\n",
        "\t| ",sch,"\n",
        "Summary:","\n",
        "\t| Size of process","\t| N  = ",format(x$N,digits=digits),".","\n",
        "\t| Number of simulation","\t| M  = ",format(x$M,digits=digits),".","\n",
        "\t| Initial values","\t| (x0,y0,z0) = ","(",format(x$x0,digits=digits),",",format(x$y0,digits=digits),",",format(x$z0,digits=digits),")",".","\n",
        "\t| Time of process","\t| t in [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t| Dt = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.snssde3d <- function(x,...) .plot.snssde3d(x,...)

lines.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z}    
    for (i in 1:x$M){
    lines(time(x),X[,i],...)
    lines(time(x),Y[,i],...)
    lines(time(x),Z[,i],...)}
}

points.snssde3d <- function(x,...)
                 {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z} 
    for (i in 1:x$M){
    points(time(x),X[,i],...)
    points(time(x),Y[,i],...)
    points(time(x),Z[,i],...)}
}


plot3D.snssde3d <- function(x,display = c("persp","rgl"),...)
                 {
	class(x) <- "snssde3d"		 
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
                  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
                  Z = x$Z} 
    plot3D(X[,1],Y[,1],Z[,1],display,...)
}

##
## summary

summary.snssde3d  <- function(object, ...)
           {   
    class(object) <- "snssde3d"
    cat("\n  Monte-Carlo Statistics for (X(t),Y(t),Z(t)) at final time T = ",object$T,"\n\n",
        sep="")
    if (object$M == 1 ){
    x <- object$X[which(time(object)==object$T)]
    y <- object$Y[which(time(object)==object$T)]
    z <- object$Z[which(time(object)==object$T)]}else{
    x <- object$X[which(time(object)==object$T),]
    y <- object$Y[which(time(object)==object$T),]
    z <- object$Z[which(time(object)==object$T),]}
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
    # rownames(res) <- paste(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          # "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  # "Bound conf Sup (95%)"),sep="")
    # names(res) <- paste(c("X","Y","Z"),sep="")
	dimnames(res) <- list(c("Mean","Variance","Median","First quartile","Third quartile","Skewness","Kurtosis","Moment of order 2",
	                          "Moment of order 3","Moment of order 4","Moment of order 5","Bound conf Inf (95%)",
							  "Bound conf Sup (95%)"),c("X","Y","Z"))
    print(res, quote = FALSE, right = TRUE,...)
    invisible(object)
}


time.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    as.vector(time(x$X))
}

mean.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    return(list(X=rowMeans(X,na.rm = TRUE,...),Y=rowMeans(Y,na.rm = TRUE,...),Z=rowMeans(Z,na.rm = TRUE,...)))
}

skewness.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    Skewx <- data.frame(sapply(1:(x$N+1),function(i) skewness(X[i,]) ))
    rownames(Skewx) <- paste("X(t=",time(x),")",sep="")
    names(Skewx) <- paste(c("Skewness"),sep="")
    Skewy <- data.frame(sapply(1:(x$N+1),function(i) skewness(Y[i,]) ))
    rownames(Skewy) <- paste("Y(t=",time(x),")",sep="")
    names(Skewy) <- paste(c("Skewness"),sep="")
    Skewz <- data.frame(sapply(1:(x$N+1),function(i) skewness(Z[i,]) ))
    rownames(Skewz) <- paste("Z(t=",time(x),")",sep="")
    names(Skewz) <- paste(c("Skewness"),sep="")	
    return(list(X=Skewx,Y=Skewy,Z=Skewz))
}

kurtosis.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    kurtx <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(X[i,]) ))
    rownames(kurtx) <- paste("X(t=",time(x),")",sep="")
    names(kurtx) <- paste(c("Kurtosis"),sep="")
    kurty <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Y[i,]) ))
    rownames(kurty) <- paste("Y(t=",time(x),")",sep="")
    names(kurty) <- paste(c("Kurtosis"),sep="")
    kurtz <- data.frame(sapply(1:(x$N+1),function(i) kurtosis(Z[i,]) ))
    rownames(kurtz) <- paste("Z(t=",time(x),")",sep="")
    names(kurtz) <- paste(c("Kurtosis"),sep="")
    return(list(X=kurtx,Y=kurty,Z=kurtz))
}

median.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    Medx <- data.frame(sapply(1:(x$N+1),function(i) median(X[i,],na.rm = TRUE) ))
    rownames(Medx) <- paste("X(t=",time(x),")",sep="")
    names(Medx) <- paste(c("Median"),sep="")
    Medy <- data.frame(sapply(1:(x$N+1),function(i) median(Y[i,],na.rm = TRUE) ))
    rownames(Medy) <- paste("Y(t=",time(x),")",sep="")
    names(Medy) <- paste(c("Median"),sep="")
    Medz <- data.frame(sapply(1:(x$N+1),function(i) median(Z[i,],na.rm = TRUE) ))
    rownames(Medz) <- paste("Z(t=",time(x),")",sep="")
    names(Medz) <- paste(c("Median"),sep="")
    return(list(X=Medx,Y=Medy,Z=Medz))
}

quantile.snssde3d <- function(x,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    Qunx <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(X[i,],na.rm = TRUE,...) ))))
    rownames(Qunx) <- paste("X(t=",time(x),")",sep="")
    Quny <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Y[i,],na.rm = TRUE,...) ))))
    rownames(Quny) <- paste("Y(t=",time(x),")",sep="")
    Qunz <- t(data.frame(do.call("cbind",lapply(1:(x$N+1),function(i) quantile(Z[i,],na.rm = TRUE,...) ))))
    rownames(Qunz) <- paste("Z(t=",time(x),")",sep="")
    return(list(X=Qunx,Y=Quny,Z=Qunz))
}

moment.snssde3d <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    Momx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(X[i,],order=order[j],...)))))
    rownames(Momx) <- paste("X(t=",time(x),")",sep="")
    names(Momx) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Y[i,],order=order[j],...)))))
    rownames(Momy) <- paste("Y(t=",time(x),")",sep="")
    names(Momy) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    Momz <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
               sapply(1:length(order), function(j) moment(Z[i,],order=order[j],...)))))
    rownames(Momz) <- paste("Z(t=",time(x),")",sep="")
    names(Momz) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(list(X=Momx,Y=Momy,Z=Momz))
}

bconfint.snssde3d <- function(x,level = 0.95,...)
                    {
    class(x) <- "snssde3d"
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)
				  Z = matrix(x$Z,nrow=length(x$Z),ncol=1)}else{
                  X = x$X
                  Y = x$Y
				  Z = x$Z}
    confx <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(confx) <- paste("X(t=",time(x),")",sep="")
    names(confx) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    confy <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Y[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(confy) <- paste("Y(t=",time(x),")",sep="")
    names(confy) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
	confz <- data.frame(do.call("rbind",lapply(1:(x$N+1), function(i) 
                      quantile(Z[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8,na.rm = TRUE) ) ) )
    rownames(confz) <- paste("Z(t=",time(x),")",sep="")
    names(confz) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(list(X=confx,Y=confy,Z=confz))
}
