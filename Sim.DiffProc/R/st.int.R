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
##### st.int

st.int <- function(expr, ...)  UseMethod("st.int")

st.int.default <- function(expr, lower = 0, upper = 1, M = 1, subdivisions = 1000L,
                          type=c("ito","str"),...)
          {
    if (any(!is.numeric(subdivisions) || (subdivisions - floor(subdivisions) > 0) || subdivisions <= 1L)) stop(" subdivisions must be a positive integer ")
    if (any(!is.numeric(M)  || (M - floor(M) > 0) || M <= 0)) stop(" 'M' must be a positive integer ")	
    if (!is.expression(expr))  stop(" 'expr' must be expression in t and w 'expr(t,w)'")
    if (missing(type)) type <- "ito"
    if (any(!is.finite(lower) && !is.finite(upper)))    stop("a limit is missing")
    if (any(is.na(lower) && is.na(upper)))              stop("a limit is missing")
    if (any(lower < 0 || upper < 0 || upper <= lower) ) stop(" limit of integration. please use positive bound 'upper > lower >= 0' ") 
    t <- seq(lower ,upper, by=(upper-lower)/subdivisions)
    fun <- function(t,w) eval(expr)
    Ito <- function()  {
            w = c(0,cumsum(rnorm(subdivisions+1,mean=0,sd=sqrt((upper-lower)/subdivisions))))
            dw   <- diff(w)
	       St <- cumsum(sapply(1:(subdivisions+1), function(i) fun(t[i],w[i])*dw[i]))
	       St    }
    Str <- function()  {
            w = c(0,cumsum(rnorm(subdivisions+1,mean=0,sd=sqrt((upper-lower)/subdivisions))))
            dw   <- diff(w)
	       St <- cumsum(sapply(1:(subdivisions+1), function(i) 0.5*(fun(t[i],w[i])+fun(t[i+1],w[i+1]))*dw[i]))
	       St    }
    if (type=="ito"){res <- data.frame(sapply(1:M,function(i) Ito()))}
    else { res <- data.frame(sapply(1:M,function(i) Str()))}
    names(res) <- paste("X",1:M,sep="")
    X <- ts(res, start = lower, deltat = (upper-lower)/subdivisions )
    structure(list(X=X, fun=expr[[1]], type=type, subdivisions=subdivisions, M = M, 
                   Dt=(upper-lower)/subdivisions,t0=lower,T=upper),class="st.int")
}

###

print.st.int <- function(x, digits=NULL, ...)
           {
    class(x) <- "st.int"
    if(x$type=="ito"){
    cat("Ito integral:","\n",
        "\t| X(t)   = integral (f(s,w) * dw(s))","\n", 
        "\t| f(t,w) = ",deparse(x$fun),"\n",
        "Summary:","\n",
        "\t| Number of subintervals","\t = ",format(x$subdivisions,digits=digits),".","\n",
        "\t| Number of simulations","\t\t = ",format(x$M,digits=digits),".","\n",
        "\t| Limits of integration","\t\t = [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t\t = ",format(x$Dt,digits=digits),".","\n",
        sep="")}else{
    cat("Stratonovich integral:","\n",
        "\t| X(t)   = integral (f(s,w) o dw(s))","\n", 
        "\t| f(t,w) = ",deparse(x$fun),"\n",
        "Summary:","\n",
        "\t| Number of subintervals","\t = ",format(x$subdivisions,digits=digits),".","\n",
        "\t| Number of simulations","\t\t = ",format(x$M,digits=digits),".","\n",
        "\t| Limits of integration","\t\t = [",format(x$t0,digits=digits),",",format(x$T,digits=digits),"].","\n",
        "\t| Discretization","\t\t = ",format(x$Dt,digits=digits),".","\n",
        sep="")}
    invisible(x)
}

##
## Plot

plot.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
    plot(X,...)
    #if(x$type=="ito"){
    #    mtext(bquote(X[t]==integral(.(x$fun)*dW[s],.(x$t0),.(x$T))),line=0.1,cex=1.1,adj=0.5)
    #}else{
    #    mtext(bquote(X[t]==integral(.(x$fun)~o~dW[s],.(x$t0),.(x$T))),line=0.1,cex=1.1,adj=0.5)
    #}
}

lines.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
    for (i in 1:dim(X)[2]){
    lines(time(x),X[,i],...)}
}

points.st.int <- function(x,...)
                 {
    class(x) <- "st.int"
    X <- x$X
    for (i in 1:dim(X)[2]){
    points(time(x),X[,i],...)}
}

add.bconfint.st.int <- function(x,level=0.95,lty=NULL,lwd=NULL,col=NULL,cex=NULL,...)
                 {
    class(x) <- "st.int"
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

summary.st.int <- function(object, ...)
           {
	x <- object	   
    class(x) <- "st.int"
    cat("\n\tMonte-Carlo Statistics for X(t) at final time T = ",x$T,"\n\n",        
        "| Process mean","\t\t\t = ",mean(x$X[which(time(x)==x$T),]),"\n",
        "| Process variance","\t\t = ",var(x$X[which(time(x)==x$T),]),"\n",
        "| Process median","\t\t = ",median(x$X[which(time(x)==x$T),]),"\n",
        "| Process first quartile","\t = ",quantile(x$X[which(time(x)==x$T),],0.25),"\n",
        "| Process third quartile","\t = ",quantile(x$X[which(time(x)==x$T),],0.75),"\n",
        "| Process skewness","\t\t = ",skewness(x$X[which(time(x)==x$T),]),"\n",
        "| Process kurtosis","\t\t = ",kurtosis(x$X[which(time(x)==x$T),]),"\n",
        "| Process moment of order 2","\t = ",moment(x$X[which(time(x)==x$T),],order=2),"\n",
        "| Process moment of order 3","\t = ",moment(x$X[which(time(x)==x$T),],order=3),"\n",
        "| Process moment of order 4","\t = ",moment(x$X[which(time(x)==x$T),],order=4),"\n",
        "| Process moment of order 5","\t = ",moment(x$X[which(time(x)==x$T),],order=5),"\n",
        "| Bound of confidence (95%)","\t = [",bconfint(x$X[which(time(x)==x$T),])[1],",",bconfint(x$X[which(time(x)==x$T),])[2],"]","\n",
        "  for the trajectories","\n",
        sep="")
    invisible(x)
}

mean.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    rowMeans(x$X,...)
}

skewness.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    Skew <- data.frame(sapply(1:(x$subdivisions+1),function(i) skewness(x$X[i,]) ))
    rownames(Skew) <- paste("X(t=",time(x),")",sep="")
    names(Skew) <- paste(c("Skewness"),sep="")
    return(Skew)
}

kurtosis.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    kurt <- data.frame(sapply(1:(x$subdivisions+1),function(i) kurtosis(x$X[i,]) ))
    rownames(kurt) <- paste("X(t=",time(x),")",sep="")
    names(kurt) <- paste(c("Kurtosis"),sep="")
    return(kurt)
}

median.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    Med <- data.frame(sapply(1:(x$subdivisions+1),function(i) median(x$X[i,]) ))
    rownames(Med) <- paste("X(t=",time(x),")",sep="")
    names(Med) <- paste(c("Median"),sep="")
    return(Med)
}

quantile.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    Qun <- t(data.frame(do.call("cbind",lapply(1:(x$subdivisions+1),function(i) quantile(x$X[i,],...) ))))
    rownames(Qun) <- paste("X(t=",time(x),")",sep="")
    return(Qun)
}

moment.st.int <- function(x,order = 2,...)
                    {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer")
    class(x) <- "st.int"
    Mom <- data.frame(do.call("rbind",lapply(1:(x$subdivisions+1), function(i) 
               sapply(1:length(order), function(j) moment(x$X[i,],order=order[j],...)))))
    rownames(Mom) <- paste("X(t=",time(x),")",sep="")
    names(Mom) <- paste(c(rep("order = ",length(order))),c(order),sep="")
    return(Mom)
}

bconfint.st.int <- function(x,level = 0.95,...)
                    {
    class(x) <- "st.int"
    conf <- data.frame(do.call("rbind",lapply(1:(x$subdivisions+1), function(i) 
                      quantile(x$X[i,], c(0.5*(1-level), 1-0.5*(1-level)),type=8) ) ) )
    rownames(conf) <- paste("X(t=",time(x),")",sep="")
    names(conf) <- paste(c(0.5*(1-level)*100,(1-(1-level)/2)*100),c(" %"," %"),sep="")
    return(conf)
}

time.st.int <- function(x,...)
                    {
    class(x) <- "st.int"
    as.vector(time(x$X))
}
