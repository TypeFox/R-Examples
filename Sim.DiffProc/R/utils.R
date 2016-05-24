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




###
### Computes the bound of the confidence

bconfint <- function(x, ...)  UseMethod("bconfint")

bconfint.default <- function(x,level = 0.95,...)
            {
   return(quantile(x,c(0.5*(1-level), 1-0.5*(1-level)),type=7,na.rm=TRUE,...))
         }
###

add.bconfint <- function(x,...) UseMethod("add.bconfint")

add.bconfint.default <- function(x,level = 0.95,...)
            {
    lines(time(x),bconfint(x,level)[,1],...)
    lines(time(x),bconfint(x,level)[,2],...)
         }

###
### Computes the moment

moment <- function(x, ...)  UseMethod("moment")

moment.default <- function(x, order = 2,...)
            {
    if (any(!is.numeric(order)  || (order - floor(order) > 0) || order < 1)) stop(" 'order' must be a positive integer ")
    ifelse(order == 1,return(0),return(mean((x - mean(x,na.rm = TRUE))^order,na.rm = TRUE)))
         }

###
### Computes the sample skewness

skewness <- function(x,...) UseMethod("skewness")

skewness.default <- function(x,...)
                   {
     return(mean((x-mean(x,na.rm = TRUE))^3,na.rm = TRUE)/sd(x,na.rm = TRUE)^3)
}

###
### Computes the sample kurtosis

kurtosis <- function(x,...) UseMethod("kurtosis")

kurtosis.default <- function(x,...)
                   {
     return(mean((x-mean(x,na.rm = TRUE))^4,na.rm = TRUE)/sd(x,na.rm = TRUE)^4)
}

###
### Plot 2D and 3D

plot2d   <- function(x,...) UseMethod("plot2d")
lines2d  <- function(x,...) UseMethod("lines2d")
points2d <- function(x,...) UseMethod("points2d")

plot2d.default <- function(x,...)
        {
    class(x) <- "plot2d"
    plot(x,...)
}

lines2d.default <- function(x,...)
        {
    class(x) <- "lines2d"
    lines(x,...)
}

points2d.default <- function(x,...)
        {
    class(x) <- "points2d"
    points(x,...)
}


plot3D   <- function(x, ...)  UseMethod("plot3D")

plot3DD <- function(X,Y,Z,display = c("persp","rgl"),col=NULL,lwd=NULL,pch=NULL,
                            type = NULL,cex=NULL,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
                            zlab=NULL,grid=NULL,angle=NULL,...)
                 {
    display <- match.arg(display)
    # if ((display == "rgl") && !(require(rgl)) ) 
                 # stop("The 'rgl' package is not available.")        
    # else if ((display == "persp") && !(require(scatterplot3d)) ) 
                 # stop("The 'scatterplot3d' package is not available.")

    if (is.null(lwd))  {lwd = 1}
    if (is.null(col))  {col = 2}
    if (is.null(pch))  {pch = 16}
    if (is.null(cex))  {cex = 0.6}
    if (is.null(type)) {type = "l"}
    if (is.null(main)) {main = ""}
    if (is.null(sub))  {sub = ""}
    if (is.null(xlab)) {xlab = expression(X[t])}
    if (is.null(ylab)) {ylab = expression(Y[t])}
    if (is.null(zlab)) {zlab = expression(Z[t])} 
    if (is.null(grid)) {grid = TRUE}   
    if (is.null(angle)){angle = 140} 
    if (display=="persp"){
         scatterplot3d(X,Y,Z,angle =angle,color=col,lwd=lwd,type=type,pch=pch,
             main = main, sub = sub,xlab = xlab, ylab = ylab, zlab = zlab,
             grid = grid,cex.symbols=cex,...)
    }else{
         plot3d(X,Y,Z,col=col,lwd=lwd,type=type,pch=pch,main = main, sub = sub,
             xlab = xlab, ylab = ylab, zlab = zlab,size=cex,...)
         }
}

plot3D.default <- function(x,display = c("persp", "rgl"),...) plot3DD(x,display,...)

####
#### plot for calss snssde2d

.plot.snssde2d <- function(x,union = TRUE,legend=TRUE,pos=1,main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde2d"
    if (is.null(col)){col=c(1,2)}
    if (is.null(lty)){lty=c(1,1)}
    if (is.null(lwd)){lwd=c(1,1)}
    if (is.null(cex)){cex=0.75}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(text.col)){text.col=c(1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    if (union){
    plot(ts.union(x$X,x$Y),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (x$M == 1){
    lines(time(x),x$X,col=col[1],lty=lty[1],lwd=lwd[1],...)
    lines(time(x),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],...)}else{
    for (i in 1:x$M){
    lines(time(x),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],...)
    lines(time(x),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],...)
         }
    }
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    }
}

.plot2d.snssde2d <- function(x,type=NULL,xlab=NULL,ylab=NULL,...)
                 {
    class(x) <- "snssde2d"
    if (is.null(type)){type="l"}
    if (is.null(ylab)){ylab=expression(Y[t])}
    if (is.null(xlab)){xlab=expression(X[t])}
    if (x$M == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)
                  Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{
                  X = x$X
                  Y = x$Y} 
    plot2d(X[,1],Y[,1],type=type,ylab=ylab,xlab=xlab,...)
    for(i in 3:4) axis(i)
}

####
#### plot for calss snssde3d

.plot.snssde3d <- function(x,union = TRUE,legend=TRUE,pos=1,main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde3d"
    if (is.null(col)){col=c(1,2,3)}
    if (is.null(lty)){lty=c(1,1,1)}
    if (is.null(lwd)){lwd=c(1,1,1)}
    if (is.null(cex)){cex=0.75}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(text.col)){text.col=c(1,1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    if (union){
    plot(ts.union(x$X,x$Y,x$Z),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (x$M == 1){
    lines(time(x),x$X,col=col[1],lty=lty[1],lwd=lwd[1],...)
    lines(time(x),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],...)
    lines(time(x),x$Z,col=col[3],lty=lty[3],lwd=lwd[3],...)}else{
    for (i in 1:x$M){
    lines(time(x),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],...)
    lines(time(x),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],...)
    lines(time(x),x$Z[,i],col=col[3],lty=lty[3],lwd=lwd[3],...)
         }
    }
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t]),expression(Z[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    dev.new()
    plot(x$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,...)
    }
}

####
#### plot for calss bridgesde2d 

.plot.bridgesde2d  <- function(x,union = TRUE,legend=TRUE,pos=1,main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "bridgesde2d"
    if (is.null(col)){col=c(1,2)}
    if (is.null(lty)){lty=c(1,1)}
    if (is.null(lwd)){lwd=c(1,1)}
    if (is.null(cex)){cex=0.75}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(text.col)){text.col=c(1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    if (union){
    plot(ts.union(x$X,x$Y),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (length(which(!is.na(x$Cx))) == 1 ){lines(as.vector(time(x$X)),x$X,col=col[1],lty=lty[1],lwd=lwd[1],...)}else{
	for (i in 1:length(which(!is.na(x$Cx)))) {lines(as.vector(time(x$X)),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],...)}}
    if (length(which(!is.na(x$Cy))) == 1 ){lines(as.vector(time(x$Y)),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],...)}else{
	for (i in 1:length(which(!is.na(x$Cy)))) {lines(as.vector(time(x$Y)),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],...)}}
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    }
}

.plot2d.bridgesde2d <- function(x,type=NULL,xlab=NULL,ylab=NULL,...)
                 {
    class(x) <- "bridgesde2d"
    if (is.null(type)){type="l"}
    if (is.null(ylab)){ylab=expression(Y[t])}
    if (is.null(xlab)){xlab=expression(X[t])}
    if (length(which(!is.na(x$Cx))) == 1){X = matrix(x$X,nrow=length(x$X),ncol=1)}else{X = x$X}  
    if (length(which(!is.na(x$Cy))) == 1){Y = matrix(x$Y,nrow=length(x$Y),ncol=1)}else{Y = x$Y}	
    plot2d(X[,1],Y[,1],type=type,ylab=ylab,xlab=xlab,...)
    for(i in 3:4) axis(i)
}

####
#### plot for calss bridgesde3d
 
.plot.bridgesde3d <- function(x,union = TRUE,legend=TRUE,pos=1,main=NULL,col=NULL,lty=NULL,lwd=NULL,cex=NULL,
                          las=NULL,text.col=NULL,...) 
              {
    class(x) <- "snssde3d"
    if (is.null(col)){col=c(1,2,3)}
    if (is.null(lty)){lty=c(1,1,1)}
    if (is.null(lwd)){lwd=c(1,1,1)}
    if (is.null(cex)){cex=0.75}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (is.null(text.col)){text.col=c(1,1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    if (union){
    plot(ts.union(x$X,x$Y,x$Z),plot.type="single",ylab="",type="n",las=las,main=main,...)
    if (length(which(!is.na(x$Cx))) == 1 ){lines(as.vector(time(x$X)),x$X,col=col[1],lty=lty[1],lwd=lwd[1],...)}else{
	for (i in 1:length(which(!is.na(x$Cx)))) {lines(as.vector(time(x$X)),x$X[,i],col=col[1],lty=lty[1],lwd=lwd[1],...)}}
    if (length(which(!is.na(x$Cy))) == 1 ){lines(as.vector(time(x$Y)),x$Y,col=col[2],lty=lty[2],lwd=lwd[2],...)}else{
	for (i in 1:length(which(!is.na(x$Cy)))) {lines(as.vector(time(x$Y)),x$Y[,i],col=col[2],lty=lty[2],lwd=lwd[2],...)}}
    if (length(which(!is.na(x$Cz))) == 1 ){lines(as.vector(time(x$Z)),x$Z,col=col[3],lty=lty[3],lwd=lwd[3],...)}else{
	for (i in 1:length(which(!is.na(x$Cz)))) {lines(as.vector(time(x$Z)),x$Z[,i],col=col[3],lty=lty[3],lwd=lwd[3],...)}}
    if (legend){	
    legend(pos,c(expression(X[t]),expression(Y[t]),expression(Z[t])),inset = .01,col=col,lty=lty,lwd=lwd,cex=cex,text.col=text.col)}
    }else{
    plot(x$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    dev.new()
    plot(x$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    dev.new()
    plot(x$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,...)
    }
}


####
#### plot for calss fptsde1d

.plot.fptsde1d <- function(x,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,type=NULL,las=NULL,ylab=NULL,...)
                 {
    class(x) <- "fptsde1d"
    Bn  <- function(t)  eval(x$boundary)+0*t
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(lty)){lty=1}
    if (is.null(las)){las=1}
    if (is.null(type)){type="l"}
	if (is.null(ylab)){ylab=expression(X[t])}
    if (is.null(text.col)){text.col=c(1,1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    plot(x$SDE,plot.type="single",col=col,lwd=lwd,lty=lty,type=type,las=las,ylab=ylab,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=2,lty=2)
    if (length(which(is.na(x$fpt)==TRUE)) < x$SDE$M){
    points(x$fpt,Bn(x$fpt),col = 4, pch = 16,cex=0.7)
    }
    ##if (x$SDE$M > 1){MM = max(x$SDE$X[x$SDE$N +1,],na.rm = TRUE)}else{MM = x$SDE$X[x$SDE$N +1]}
    ##if (MM >= x$SDE$x0){pos = "topleft"}else{pos = "topright"} 
    if (legend) {
    Dr <- gsub(pattern = 'x', replacement = 'X[t]', x = as.expression(x$SDE$drift), ignore.case = F,fixed = T)
    DD <- gsub(pattern = 'x', replacement = 'X[t]', x = as.expression(x$SDE$diffusion), ignore.case = F,fixed = T)
    if (x$SDE$type == "ito"){
    A = as.expression(bquote(SDE: dX[t] == .(parse(text=Dr,srcfile=Dr)[[1]]) *~dt + .(parse(text=DD,srcfile=DD)[[1]])~dW[t]))}else{
    A = as.expression(bquote(SDE: dX[t] == .(parse(text=Dr,srcfile=Dr)[[1]]) *~dt + .(parse(text=DD,srcfile=DD)[[1]])~o~dW[t]))}
    B = as.expression(bquote(Boundary: S[t]==.(x$boundary)))
    if (x$SDE$x0 > Bn(x$SDE$t0)){
    FPT = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] <= S[t]),"}")))}else{
    FPT = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] >= S[t]),"}")))}
    if (length(which(is.na(x$fpt))) == x$SDE$M ){
    legend(pos,legend=c(A,B), inset = .01,lty = c(lty, 2),pch=c(NA,NA),
           lwd=c(lwd,1),col=c(col,2),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A,B,FPT), inset = .01,lty = c(lty, 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd,1,NA),col=c(col,2,4),text.col=text.col,cex=cex)
        }
    }
}

####
#### plot for calss fptsde2d

.plot.fptsde2d <- function(x,union = TRUE,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,las=NULL,main=NULL,...)
                 {
    class(x) <- "fptsde2d"
    Bn  <- function(t)  eval(x$boundary)+0*t
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1:2}
    if (is.null(lwd)){lwd=c(1,1)}
    if (is.null(lty)){lty=c(1,1)}
    if (is.null(las)){las=1}
    if (is.null(text.col)){text.col=c(1,1,1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
	A1 <-expression(X[t])
	A2 <-expression(Y[t])
    B = as.expression(bquote(S[t]==.(x$boundary)))
    if (x$SDE$x0 > Bn(x$SDE$t0)){
    FPTx = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] <= S[t]),"}")))}else{
    FPTx = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] >= S[t]),"}")))}
    if (x$SDE$y0 > Bn(x$SDE$t0)){
    FPTy = as.expression(bquote(Tau[(list(Y[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Y[t] <= S[t]),"}")))}else{
    FPTy = as.expression(bquote(Tau[(list(Y[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Y[t] >= S[t]),"}")))}
    if (union){
    plot(x$SDE,legend=FALSE,main=main,col=col,lty=lty,lwd=lwd,cex=cex,las=las,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=3,lty=2)
    if (length(which(is.na(x$fptx)==TRUE)) < x$SDE$M){
    points(x$fptx,Bn(x$fptx),col = 4, pch = 16,cex=0.7)
    }
    if (length(which(is.na(x$fpty)==TRUE)) < x$SDE$M){
    points(x$fpty,Bn(x$fpty),col = 5, pch = 16,cex=0.7)
    }
    if (legend) {
	Drx <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = as.expression(x$SDE$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDx <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = as.expression(x$SDE$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	Dry <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = as.expression(x$SDE$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDy <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = as.expression(x$SDE$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    ##if (x$SDE$type == "ito"){
    ##A1 = as.expression(bquote(dX[t] == .(parse(text=Drx,srcfile=Drx)[[1]]) *~dt + .(parse(text=DDx,srcfile=DDx)[[1]])~dW[t]^1))
    ##A2 = as.expression(bquote(dY[t] == .(parse(text=Dry,srcfile=Dry)[[1]]) *~dt + .(parse(text=DDy,srcfile=DDy)[[1]])~dW[t]^2))}else{
    ##A1 = as.expression(bquote(dX[t] == .(parse(text=Drx,srcfile=Drx)[[1]]) *~dt + .(parse(text=DDx,srcfile=DDx)[[1]])~o~dW[t]^1))
    ##A2 = as.expression(bquote(dY[t] == .(parse(text=Dry,srcfile=Dry)[[1]]) *~dt + .(parse(text=DDy,srcfile=DDy)[[1]])~o~dW[t]^2))}
    if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) != x$SDE$M){
    legend(pos,legend=c(A1,A2,B,FPTy), inset = .01,lty = c(lty[1],lty[2], 2,NA),pch=c(NA,NA,NA,16),
           lwd=c(lwd[1],lwd[2],1,NA),col=c(col[1],col[2],3,5),text.col=text.col,cex=cex)}
	else if (length(which(is.na(x$fptx))) != x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M ){
    legend(pos,legend=c(A1,A2,B,FPTx), inset = .01,lty = c(lty[1],lty[2], 2,NA),pch=c(NA,NA,NA,16),
           lwd=c(lwd[1],lwd[2],1,NA),col=c(col[1],col[2],3,4),text.col=text.col,cex=cex)}		   
	else if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M ){
    legend(pos,legend=c(A1,A2,B), inset = .01,lty = c(lty[1], lty[2], NA),pch=c(NA,NA,NA),
           lwd=c(lwd[1],lwd[2],1),col=c(col[1],col[2],3),text.col=text.col,cex=cex)}
    else {legend(pos,legend=c(A1,A2,B,FPTx,FPTy), inset = .01,lty = c(lty[1],lty[2], 2, NA, NA),pch=c(NA,NA,NA,16,16),
                lwd=c(lwd[1],lwd[2],1,NA,NA),col=c(col[1],col[2],3,4,5),text.col=text.col,cex=cex)}		   
	 }
	}else{
    plot(x$SDE$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
	lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=3,lty=2)
    if (length(which(is.na(x$fptx)==TRUE)) < x$SDE$M){
    points(x$fptx,Bn(x$fptx),col = 4, pch = 16,cex=0.7)
    }
	if (legend){
	if (length(which(is.na(x$fptx))) == x$SDE$M ){
    legend(pos,legend=c(A1,B), inset = .01,lty = c(lty[1], 2),pch=c(NA,NA),
           lwd=c(lwd[1],1),col=c(col[1],3),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A1,B,FPTx), inset = .01,lty = c(lty[1], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[1],1,NA),col=c(col[1],3,4),text.col=text.col,cex=cex)
            }
		}
    dev.new()
    plot(x$SDE$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=3,lty=2)
    if (length(which(is.na(x$fpty)==TRUE)) < x$SDE$M){
    points(x$fpty,Bn(x$fpty),col = 5, pch = 16,cex=0.7)
    }
	if (legend){
	if (length(which(is.na(x$fpty))) == x$SDE$M ){
    legend(pos,legend=c(A2,B), inset = .01,lty = c(lty[2], 2),pch=c(NA,NA),
           lwd=c(lwd[2],1),col=c(col[2],3),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A2,B,FPTy), inset = .01,lty = c(lty[2], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[2],1,NA),col=c(col[2],3,5),text.col=text.col,cex=cex)
            }
		}
	}
}	

####
#### plot for calss fptsde3d

.plot.fptsde3d <- function(x,union = TRUE,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,las=NULL,main=NULL,...)
                 {
    class(x) <- "fptsde3d"
    Bn  <- function(t)  eval(x$boundary)+0*t
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1:3}
    if (is.null(lwd)){lwd=c(1,1,1)}
    if (is.null(lty)){lty=c(1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(text.col)){text.col=c(1,1,1,1,1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
	A1 <-expression(X[t])
	A2 <-expression(Y[t])
    A3 <-expression(Z[t])
    B = as.expression(bquote(S[t]==.(x$boundary)))
    if (x$SDE$x0 > Bn(x$SDE$t0)){
    FPTx = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] <= S[t]),"}")))}else{
    FPTx = as.expression(bquote(Tau[(list(X[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),X[t] >= S[t]),"}")))}
    if (x$SDE$y0 > Bn(x$SDE$t0)){
    FPTy = as.expression(bquote(Tau[(list(Y[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Y[t] <= S[t]),"}")))}else{
    FPTy = as.expression(bquote(Tau[(list(Y[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Y[t] >= S[t]),"}")))}
    if (x$SDE$z0 > Bn(x$SDE$t0)){
    FPTz = as.expression(bquote(Tau[(list(Z[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Z[t] <= S[t]),"}")))}else{
    FPTz = as.expression(bquote(Tau[(list(Z[t],S[t]))] == inf*group("{",list(t >= .(x$SDE$t0),Z[t] >= S[t]),"}")))}	
    if (union){
    plot(x$SDE,legend=FALSE,main=main,col=col,lty=lty,lwd=lwd,cex=cex,las=las,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=4,lty=2)
    if (length(which(is.na(x$fptx)==TRUE)) < x$SDE$M){
    points(x$fptx,Bn(x$fptx),col = 5, pch = 16,cex=0.7)
    }
    if (length(which(is.na(x$fpty)==TRUE)) < x$SDE$M){
    points(x$fpty,Bn(x$fpty),col = 6, pch = 16,cex=0.7)
    }
    if (length(which(is.na(x$fptz)==TRUE)) < x$SDE$M){
    points(x$fptz,Bn(x$fptz),col = 7, pch = 16,cex=0.7)
    }
    if (legend) {
    Drx <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$driftx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDx <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$diffx), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    Dry <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$drifty), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDy <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$diffy), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	Drz <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$driftz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
	DDz <- gsub(pattern = 'x', replacement = 'X[t]', x = gsub(pattern = 'y', replacement = 'Y[t]', x = gsub(pattern = 'z', replacement = 'Z[t]', x = as.expression(x$diffz), ignore.case = F,fixed = T), ignore.case = F,fixed = T), ignore.case = F,fixed = T)
    if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) != x$SDE$M & length(which(is.na(x$fptz))) != x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTy,FPTz), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA,NA),pch=c(NA,NA,NA,NA,16,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA,NA),col=c(col[1],col[2],col[3],4,6,7),text.col=text.col,cex=cex)}
	else if (length(which(is.na(x$fptx))) != x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M & length(which(is.na(x$fptz))) != x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTx,FPTz), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA,NA),pch=c(NA,NA,NA,NA,16,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA,NA),col=c(col[1],col[2],col[3],4,5,7),text.col=text.col,cex=cex)}	
	else if (length(which(is.na(x$fptx))) != x$SDE$M & length(which(is.na(x$fpty))) != x$SDE$M & length(which(is.na(x$fptz))) == x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTx,FPTy), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA,NA),pch=c(NA,NA,NA,NA,16,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA,NA),col=c(col[1],col[2],col[3],4,5,6),text.col=text.col,cex=cex)}	
	else if (length(which(is.na(x$fptx))) != x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M & length(which(is.na(x$fptz))) == x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTx), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA),pch=c(NA,NA,NA,NA,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA),col=c(col[1],col[2],col[3],4,5),text.col=text.col,cex=cex)}	
	else if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) != x$SDE$M & length(which(is.na(x$fptz))) == x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTy), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA),pch=c(NA,NA,NA,NA,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA),col=c(col[1],col[2],col[3],4,6),text.col=text.col,cex=cex)}		
	else if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M & length(which(is.na(x$fptz))) != x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B,FPTz), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA),pch=c(NA,NA,NA,NA,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA),col=c(col[1],col[2],col[3],4,7),text.col=text.col,cex=cex)}	
	else if (length(which(is.na(x$fptx))) == x$SDE$M & length(which(is.na(x$fpty))) == x$SDE$M & length(which(is.na(x$fptz))) == x$SDE$M){
    legend(pos,legend=c(A1,A2,A3,B), inset = .01,lty = c(lty[1],lty[2],lty[3], 2),pch=c(NA,NA,NA,NA),
           lwd=c(lwd[1],lwd[2],lwd[3],1),col=c(col[1],col[2],col[3],4),text.col=text.col,cex=cex)}	
	else {
    legend(pos,legend=c(A1,A2,A3,B,FPTx,FPTy,FPTz), inset = .01,lty = c(lty[1],lty[2],lty[3], 2,NA,NA,NA),pch=c(NA,NA,NA,NA,16,16,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA,NA,NA),col=c(col[1],col[2],col[3],4,5,6,7),text.col=text.col,cex=cex)}		   	   
	 }
	}else{
    plot(x$SDE$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
	lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=4,lty=2)
    if (length(which(is.na(x$fptx)==TRUE)) < x$SDE$M){
    points(x$fptx,Bn(x$fptx),col = 5, pch = 16,cex=0.7)
    }
	if (legend){
	if (length(which(is.na(x$fptx))) == x$SDE$M ){
    legend(pos,legend=c(A1,B), inset = .01,lty = c(lty[1], 2),pch=c(NA,NA),
           lwd=c(lwd[1],1),col=c(col[1],4),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A1,B,FPTx), inset = .01,lty = c(lty[1], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[1],1,NA),col=c(col[1],4,5),text.col=text.col,cex=cex)
            }
		}
    dev.new()
    plot(x$SDE$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=4,lty=2)
    if (length(which(is.na(x$fpty)==TRUE)) < x$SDE$M){
    points(x$fpty,Bn(x$fpty),col = 6, pch = 16,cex=0.7)
    }
	if (legend){
	if (length(which(is.na(x$fpty))) == x$SDE$M ){
    legend(pos,legend=c(A2,B), inset = .01,lty = c(lty[2], 2),pch=c(NA,NA),
           lwd=c(lwd[2],1),col=c(col[2],4),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A2,B,FPTy), inset = .01,lty = c(lty[2], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[2],1,NA),col=c(col[2],4,6),text.col=text.col,cex=cex)
            }
		}
    dev.new()
    plot(x$SDE$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,...)
    lines(seq(x$SDE$t0,x$SDE$T,length.out=1000),Bn(seq(x$SDE$t0,x$SDE$T,length.out=1000)),col=4,lty=2)
    if (length(which(is.na(x$fptz)==TRUE)) < x$SDE$M){
    points(x$fptz,Bn(x$fptz),col = 7, pch = 16,cex=0.7)
    }
	if (legend){
	if (length(which(is.na(x$fptz))) == x$SDE$M ){
    legend(pos,legend=c(A3,B), inset = .01,lty = c(lty[3], 2),pch=c(NA,NA),
           lwd=c(lwd[3],1),col=c(col[3],4),text.col=text.col,cex=cex)}else{
    legend(pos,legend=c(A3,B,FPTz), inset = .01,lty = c(lty[3], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[3],1,NA),col=c(col[3],4,7),text.col=text.col,cex=cex)
            }
		}	
	}
}	

####
#### plot for calss rsde1d

.plot.rsde1d <- function(x,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,type=NULL,las=NULL,ylab=NULL,...)
                 {
    class(x) <- "rsde1d"
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1}
    if (is.null(lwd)){lwd=1}
    if (is.null(lty)){lty=1}
    if (is.null(las)){las=1}
    if (is.null(type)){type="l"}
	if (is.null(ylab)){ylab=expression(X[t])}
    if (is.null(text.col)){text.col=c(1,1,1)}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    plot(x$SDE,plot.type="single",col=col,lwd=lwd,lty=lty,type=type,las=las,ylab=ylab,...)
    points(rep(x$tau,x$SDE$M),x$x,col = 3, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=2, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$x,na.rm = TRUE)),lty=2,col=2)
    lines(c(x$tau,x$tau),c(0,min(x$x,na.rm = TRUE)),lty=2,col=2)
    if (legend) {
    A = expression(X[t])
    B = bquote(tau ==.(x$tau))
    X = bquote(x[tau] == group("{",list(t >= .(x$SDE$t0),X[t] == X[.(x$tau)]),"}"))
    legend(pos,legend=c(A,B,X), inset = .01,lty = c(lty, 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd,1,NA),col=c(col,2,3),text.col=text.col,cex=cex)
        }
}

####
#### plot for calss rsde2d

.plot.rsde2d <- function(x,union = TRUE,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,las=NULL,main=NULL,...)
                 {
    class(x) <- "rsde2d"
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1:2}
    if (is.null(lwd)){lwd=c(1,1)}
    if (is.null(lty)){lty=c(1,1)}
    if (is.null(las)){las=1}
    if (is.null(text.col)){text.col=c(1,1,1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    A1 = expression(X[t])
	A2 = expression(Y[t])
    B = bquote(tau ==.(x$tau))
    X = bquote(x[tau] == group("{",list(t >= .(x$SDE$t0),X[t] == X[.(x$tau)]),"}"))	
    Y = bquote(y[tau] == group("{",list(t >= .(x$SDE$t0),Y[t] == Y[.(x$tau)]),"}"))	
    if (union){
    plot(x$SDE,legend=FALSE,main=main,col=col,lty=lty,lwd=lwd,cex=cex,las=las,...)
    points(rep(x$tau,x$SDE$M),x$x,col = 4, pch = 16,cex=0.8)
	points(rep(x$tau,x$SDE$M),x$y,col = 5, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=3, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(max(x$x,na.rm = TRUE),max(x$y,na.rm = TRUE))),lty=2,col=3)
    lines(c(x$tau,x$tau),c(0,min(min(x$x,na.rm = TRUE),min(x$y,na.rm = TRUE))),lty=2,col=3)
    if (legend) {
    legend(pos,legend=c(A1,A2,B,X,Y), inset = .01,lty = c(lty[1],lty[2], 2, NA,NA),pch=c(NA,NA,NA,16,16),
           lwd=c(lwd[1],lwd[2],1,NA,NA),col=c(col[1],col[2],3,4,5),text.col=text.col,cex=cex)	   
	 } 
	}else{
    plot(x$SDE$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    points(rep(x$tau,x$SDE$M),x$x,col = 4, pch = 16,cex=0.8)
	Axis(at = x$tau, side=1,col=3, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$x,na.rm = TRUE)),lty=2,col=3)
    lines(c(x$tau,x$tau),c(0,min(x$x,na.rm = TRUE)),lty=2,col=3)
	if (legend){
    legend(pos,legend=c(A1,B,X), inset = .01,lty = c(lty[1], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[1],1,NA),col=c(col[1],3,4),text.col=text.col,cex=cex)
		}
    dev.new()
    plot(x$SDE$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
	points(rep(x$tau,x$SDE$M),x$y,col = 5, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=3, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$y,na.rm = TRUE)),lty=2,col=3)
    lines(c(x$tau,x$tau),c(0,min(x$y,na.rm = TRUE)),lty=2,col=3)
	if (legend){
    legend(pos,legend=c(A2,B,Y), inset = .01,lty = c(lty[2], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[2],1,NA),col=c(col[2],3,5),text.col=text.col,cex=cex)
		}
	}
}	

####
#### plot for calss rsde3d

.plot.rsde3d <- function(x,union = TRUE,legend=TRUE,pos=2,cex=NULL,col=NULL,lwd=NULL,
                          text.col=NULL,lty=NULL,las=NULL,main=NULL,...)
                 {
    class(x) <- "rsde3d"
    if (is.null(cex)){cex=0.72}
    if (is.null(col)){col=1:3}
    if (is.null(lwd)){lwd=c(1,1,1)}
    if (is.null(lty)){lty=c(1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(text.col)){text.col=c(1,1,1,1,1,1,1)}
    if (is.null(las)){las=1}
    if (is.null(main)){main=""}
    if (pos==1){pos = "top"}
    else if (pos==2){pos = "topright"}
    else if (pos==3){pos = "topleft"}
    else if (pos==4){pos = "center"}
    else if (pos==5){pos = "right"}
    else if (pos==6){pos = "left"}
    else if (pos==7){pos = "bottom"}
    else if (pos==8){pos = "bottomright"}
    else if (pos==9){pos = "bottomleft"}
    A1 = expression(X[t])
	A2 = expression(Y[t])
	A3 = expression(Z[t])
    B = bquote(tau ==.(x$tau))
    X = bquote(x[tau] == group("{",list(t >= .(x$SDE$t0),X[t] == X[.(x$tau)]),"}"))	
    Y = bquote(y[tau] == group("{",list(t >= .(x$SDE$t0),Y[t] == Y[.(x$tau)]),"}"))	
    Z = bquote(z[tau] == group("{",list(t >= .(x$SDE$t0),Z[t] == Z[.(x$tau)]),"}"))	
    if (union){
    plot(x$SDE,legend=FALSE,main=main,col=col,lty=lty,lwd=lwd,cex=cex,las=las,...)
    points(rep(x$tau,x$SDE$M),x$x,col = 5, pch = 16,cex=0.8)
	points(rep(x$tau,x$SDE$M),x$y,col = 6, pch = 16,cex=0.8)
	points(rep(x$tau,x$SDE$M),x$z,col = 7, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=3, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(max(x$x,na.rm = TRUE),max(x$y,na.rm = TRUE))),lty=2,col=4)
    lines(c(x$tau,x$tau),c(0,min(min(x$x,na.rm = TRUE),min(x$y,na.rm = TRUE))),lty=2,col=4)
    lines(c(x$tau,x$tau),c(0,min(min(x$z,na.rm = TRUE),min(x$z,na.rm = TRUE))),lty=2,col=4)
    if (legend) {
    legend(pos,legend=c(A1,A2,A3,B,X,Y,Z), inset = .01,lty = c(lty[1],lty[2],lty[3], 2, NA,NA,NA),pch=c(NA,NA,NA,NA,16,16,16),
           lwd=c(lwd[1],lwd[2],lwd[3],1,NA,NA,NA),col=c(col[1],col[2],col[3],4,5,6,7),text.col=text.col,cex=cex)	   
	 } 
	}else{
    plot(x$SDE$X,plot.type="single",ylab=expression(X[t]),col=col[1],lty=lty[1],lwd=lwd[1],las=las,main=main,...)
    points(rep(x$tau,x$SDE$M),x$x,col = 5, pch = 16,cex=0.8)
	Axis(at = x$tau, side=1,col=4, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$x,na.rm = TRUE)),lty=2,col=4)
    lines(c(x$tau,x$tau),c(0,min(x$x,na.rm = TRUE)),lty=2,col=4)
	if (legend){
    legend(pos,legend=c(A1,B,X), inset = .01,lty = c(lty[1], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[1],1,NA),col=c(col[1],4,5),text.col=text.col,cex=cex)
		}
    dev.new()
    plot(x$SDE$Y,plot.type="single",ylab=expression(Y[t]),col=col[2],lty=lty[2],lwd=lwd[2],las=las,main=main,...)
	points(rep(x$tau,x$SDE$M),x$y,col = 6, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=4, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$y,na.rm = TRUE)),lty=2,col=4)
    lines(c(x$tau,x$tau),c(0,min(x$y,na.rm = TRUE)),lty=2,col=4)
	if (legend){
    legend(pos,legend=c(A2,B,Y), inset = .01,lty = c(lty[2], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[2],1,NA),col=c(col[2],4,6),text.col=text.col,cex=cex)
		}
    dev.new()
    plot(x$SDE$Z,plot.type="single",ylab=expression(Z[t]),col=col[3],lty=lty[3],lwd=lwd[3],las=las,main=main,...)
	points(rep(x$tau,x$SDE$M),x$z,col = 7, pch = 16,cex=0.8)
    Axis(at = x$tau, side=1,col=4, labels = bquote(tau))
    lines(c(x$tau,x$tau),c(0,max(x$z,na.rm = TRUE)),lty=2,col=4)
    lines(c(x$tau,x$tau),c(0,min(x$z,na.rm = TRUE)),lty=2,col=4)
	if (legend){
    legend(pos,legend=c(A3,B,Z), inset = .01,lty = c(lty[3], 2, NA),pch=c(NA,NA,16),
           lwd=c(lwd[3],1,NA),col=c(col[3],4,7),text.col=text.col,cex=cex)
		}
	}
}	


####
#### Char2expression

.Char2exp <- function(expr)
         {
y <- gsub(pattern = 'theta[1]', replacement = 'theta1', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[2]', replacement = 'theta2', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[3]', replacement = 'theta3', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[4]', replacement = 'theta4', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[5]', replacement = 'theta5', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[6]', replacement = 'theta6', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[7]', replacement = 'theta7', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[8]', replacement = 'theta8', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[9]', replacement = 'theta9', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[10]', replacement = 'theta10', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[11]', replacement = 'theta11', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[12]', replacement = 'theta12', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[13]', replacement = 'theta13', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[14]', replacement = 'theta14', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[15]', replacement = 'theta15', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[16]', replacement = 'theta16', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[17]', replacement = 'theta17', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[18]', replacement = 'theta18', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[19]', replacement = 'theta19', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta[20]', replacement = 'theta20', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}


.Exp2char <- function(expr)
         {
y <- gsub(pattern = 'theta1', replacement = 'theta[1]', x = expr, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta2', replacement = 'theta[2]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta3', replacement = 'theta[3]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta4', replacement = 'theta[4]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta5', replacement = 'theta[5]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta6', replacement = 'theta[6]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta7', replacement = 'theta[7]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta8', replacement = 'theta[8]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta9', replacement = 'theta[9]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta10', replacement = 'theta[10]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta11', replacement = 'theta[11]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta12', replacement = 'theta[12]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta13', replacement = 'theta[13]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta14', replacement = 'theta[14]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta15', replacement = 'theta[15]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta16', replacement = 'theta[16]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta17', replacement = 'theta[17]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta18', replacement = 'theta[18]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta19', replacement = 'theta[19]', x = y, ignore.case = F,fixed = T)
y <- gsub(pattern = 'theta20', replacement = 'theta[20]', x = y, ignore.case = F,fixed = T)
Y <- parse(text=y,srcfile=y)
return(Y)
}

#####
##### Simulate a Multivariate Normal Distribution

.rMnorm <- function(n , Mu, Sigma) 
         {
    d <- dim(Sigma)[1]
    Egn <- eigen(Sigma)
    Y <- matrix(rnorm(d * n), n)
    X <- drop(Mu) + Egn$vectors %*% diag(sqrt(pmax(Egn$values, 0)), d) %*% t(Y)
    return(t(X))
}

#####
##### dc.sde1d

.dcEuler <- function(x,t,x0,t0,theta,drift,diffusion,log=FALSE)
              {
	 A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     dt  <- t-t0
     lik <- dnorm(x,mean=x0+A(t0,x0,theta)*dt,sd=sqrt(dt)*S(t0,x0,theta),log=log)
     lik[is.infinite(lik)] <- NA
     lik
       }
	   
.dcOzaki <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
             {  
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
     Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
     dt  <- t-t0
     K   <- log(1+(A(t0,x0,theta)/(x0*Ax(t0,x0,theta)))*(exp(Ax(t0,x0,theta)*dt)-1))/dt
     E   <- x0 + (A(t0,x0,theta)/Ax(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1)
     V   <- S(t0,x0,theta)^2 * ((exp(2*K*dt) -1)/(2*K))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log)
     lik[is.infinite(lik)] <- NA
     lik      
}

.dcShoji <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
            {
     A   <- function(t,x,theta)  eval(drift)
     S   <- function(t,x,theta)  eval(diffusion)
	 At  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"t"))))
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))
     dt  <- t-t0
     E   <- x0 + A(t0,x0,theta)*(exp(Ax(t0,x0,theta)*dt)-1)/Ax(t0,x0,theta) + (0.5 * S(t0,x0,theta)^2 * Axx(t0,x0,theta) + At(t0,x0,theta))*(exp(Ax(t0,x0,theta)*dt)-1-Ax(t0,x0,theta)*dt)/Ax(t0,x0,theta)^2 
     V   <- S(t0,x0,theta)^2 * (exp(2*Ax(t0,x0,theta)*dt)-1)/(2*Ax(t0,x0,theta))
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik   
}

# .dcElerian <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
            # {
     # test <- .Exp2char(as.expression(D(.Char2exp(diffusion),"x")))[[1]]
     # if (test == 0) stop("The approximation is not valid, because 'deriv(diffusion,x) = 0'")
     # A   <- function(t,x,theta)  eval(drift)
     # S   <- function(t,x,theta)  eval(diffusion)
	 # Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))	 
     # dt <- t-t0
     # E <- 0.5*S(t0, x0, theta)*Sx(t0, x0, theta)*dt
     # B <- -0.5* (S(t0, x0,theta)/Sx(t0, x0,theta)) + x0 + A(t0, x0, theta)*dt - E
     # C <- 1/((S(t0, x0,theta)^2) * dt)
     # z <- (x-B)/E
     # z[z<0] <- NA
     # lik <- ( exp(-0.5*(C+z)) * cosh(sqrt(C*z)) )/( sqrt(z) *abs(E)* sqrt(2*pi) )	 
     # if(log) lik <- log(lik)
     # lik[is.infinite(lik)] <- NA
     # lik
# }

.dcKessler <- function(x, t, x0, t0, theta, drift,diffusion, log=FALSE)
           {
     A   <- function(t,x,theta)  eval(drift)
	 Ax  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(drift),"x"))))
	 Axx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(drift),"x"),"x"))))	
     S   <- function(t,x,theta)  eval(diffusion)
	 Sx  <- function(t,x,theta)  eval(.Exp2char(as.expression(D(.Char2exp(diffusion),"x"))))
	 Sxx <- function(t,x,theta)  eval(.Exp2char(as.expression(D(D(.Char2exp(diffusion),"x"),"x"))))		   
     dt <- t-t0 
     E   <- x0 + A(t0,x0,theta)*dt + 0.5*(dt)^2 * (A(t0,x0,theta)*Ax(t0,x0,theta)+0.5*(S(t0,x0,theta))^2 * Axx(t0,x0,theta))
     V   <- x0^2 + (2*A(t0,x0,theta)*x0+(S(t0,x0,theta))^2)*dt + 0.5*(dt)^2 *(2*A(t0,x0,theta)*(Ax(t0,x0,theta)*x0+A(t0,x0,theta)+S(t0,x0,theta)*Sx(t0,x0,theta))+(S(t0,x0,theta))^2 *(Axx(t0,x0,theta)*x0+2*Ax(t0,x0,theta)+(Sx(t0,x0,theta))^2 + S(t0,x0,theta)*Sxx(t0,x0,theta))) - E^2
	 V[V < 0] <- NA
     lik <- dnorm(x, mean=E, sd=sqrt(V),log=log) 
	 lik[is.infinite(lik)] <- NA
     lik 
}


