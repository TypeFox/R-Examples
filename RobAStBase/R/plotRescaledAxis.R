## helper functions for rescaling x and y axis in various diagnostic plots

.rescalefct <- function(x, fct,
         scaleX = FALSE, scaleX.fct, scaleX.inv,
         scaleY = FALSE, scaleY.fct = pnorm,
         xlim, ylim, dots){

# if scaleX rescales x, if scaleY rescales fct(x);
# to this end uses trafos scaleX.fct with inverse scale.inv
# resp. scaleY.fct; it respects xlim and  ylim (given in orig. scale)
# thins out the scaled values if necessary and accordingly modifies
# slots xaxt, yaxt, axes of dots to indicate the new axes have to be drawn
#    paradigm small letters = orig. scale, capital letters = transformed scale
# return value: list with (thinned out) x and y, X and Y and modified dots

         X <- x
         wI <- 1:length(x)
         if(scaleX){
            if(!is.null(xlim)){
                   dots$xlim <- scaleX.fct(xlim)
                   x <- x[x>=xlim[1] & x<=xlim[2]]
            }
            Xo <- X <- scaleX.fct(x)
            X <- .DistrCollapse(X, 0*X)$supp
            wI <- sapply(X, function(uu){ w<- which(uu==Xo); if(length(w)>0) w[1] else NA})
            wI <- wI[!is.na(wI)]
            x <- scaleX.inv(X)
            dots$axes <- NULL
            dots$xaxt <- "n"
         }
         Y <- y <- if(is.function(fct)) fct(x) else fct[wI,1]
         scy <- if(is.function(fct)) NA else fct[wI,2]
         if(scaleY){
            Y <- scaleY.fct(y)
            if(!is.null(ylim)) dots$ylim <- scaleY.fct(ylim)
            dots$axes <- NULL
            dots$yaxt <- "n"
            }
         return(list(x=x,y=y,X=X,Y=Y,scy=scy,dots=dots))
}

if(FALSE){
  set.seed(1); x<- sort(rnorm(10))
  res <- .rescalefct(x, fct=function(s) sin(s), xlim=c(-2,1),ylim=c(0,1),dots=list(NULL))
  res2 <- .rescalefct(x, fct=function(s) sin(s), scaleY=T, xlim=c(-2,1),ylim=c(0,1),dots=list(NULL))
  res3 <- .rescalefct(x, fct=function(s) sin(s),
            scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T,
            xlim=c(-2,1),ylim=c(0,1),dots=list(NULL))
  distroptions("DistrResolution"=0.05)
  res4 <- .rescalefct(x, fct=function(s) sin(s),
            scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T,
            xlim=c(-2,1),ylim=c(0,1),dots=list(NULL))
}

.plotRescaledAxis <- function(scaleX,scaleX.fct, scaleX.inv,
                              scaleY,scaleY.fct, scaleY.inv,
                              xlim, ylim, X, ypts = 400, n = 11,
                              x.ticks = NULL, y.ticks = NULL, withbox = TRUE){
# plots rescaled axes acc. to logicals scaleX, scaleY
# to this end uses trafos scaleX.fct with inverse scale.inv
# resp. scaleY.fct; it respects xlim and  ylim (given in orig. scale)
# return value: none
        if(scaleX){
           if(is.null(x.ticks)){
               x <- pretty(scaleX.inv(X))
               if(!is.null(xlim)) x <- pmax(x, xlim[1])
               if(!is.null(xlim)) x <- pmin(x, xlim[2])
               X <- .DistrCollapse(scaleX.fct(x),0*x)$supp
               x <- scaleX.inv(X)
               x <- x[is.finite(x)]
               x <- pretty(x,n=n)
               X <- .DistrCollapse(scaleX.fct(x),0*x)$supp
               x <- scaleX.inv(X)
               x <- x[is.finite(x)]
               x <- pretty(x,n=length(x))
               x[.isEqual01(x)&x<0.4] <- 0
               X <- scaleX.fct(x)
               xf <- prettyNum(x)
               i01 <- !.isEqual01(X)
               xf <- xf[i01]
               Xi <- X
               X <- X[i01]
               i0 <- any(!i01&Xi<0.5)
               i1 <- any(!i01&Xi>0.5)
               if(i0){ xf <- c(NA,xf); X <- c(0, X)}
               if(i1){ xf <- c(xf,NA); X <- c(X, 1)}
               axis(1,at=X,labels=xf)
               if(i0) axis(1,at=0,labels=expression(-infinity))
               if(i1) axis(1,at=1,labels=expression(infinity))
            }else{
               if(is.null(xlim)){ xlim <- c(-Inf,Inf)}else{
                  if(is.na(xlim[1])) xlim[1] <- -Inf
                  if(is.na(xlim[2])) xlim[2] <- Inf }
               x.ticks <- sort(unique(x.ticks[!is.na(x.ticks)]))
               xf <- pmin(pmax(x.ticks[is.finite(x.ticks)],xlim[1]),xlim[2])
               Xf <- scaleX.fct(xf)
               axis(1,at=Xf,labels=xf)
               if(-Inf %in% x.ticks) axis(1,at=0,labels=expression(-infinity))
               if(Inf %in% x.ticks)  axis(1,at=1,labels=expression(infinity))
            }
            box()
        }else{
            if(!is.null(x.ticks)){
               if(is.null(xlim)){ xlim <- c(-Inf,Inf)}else{
                  if(is.na(xlim[1])) xlim[1] <- -Inf
                  if(is.na(xlim[2])) xlim[2] <- Inf }
               x.ticks <- sort(unique(x.ticks[!is.na(x.ticks)]))
               xf <- pmin(pmax(x.ticks[is.finite(x.ticks)],xlim[1]),xlim[2])
               axis(1,at=xf,labels=xf)
               if(-Inf %in% x.ticks) axis(1,at=0,labels=expression(-infinity))
               if(Inf %in% x.ticks)  axis(1,at=1,labels=expression(infinity))
               box()
            }
        }
        if(scaleY){
           if(is.null(y.ticks)){
               Y0 <- if(!is.null(ylim)) max(0, scaleY.fct(ylim[1])) else 0
               Y1 <- if(!is.null(ylim)) min(1, scaleY.fct(ylim[2])) else 1
               Y <- seq(Y0,Y1, length=ypts)
               y <- pretty(scaleY.inv(Y),n=n)
               Y <- .DistrCollapse(scaleY.fct(y),0*y)$supp
               y <- scaleY.inv(Y)
               y <- y[is.finite(y)]
               y <- pretty(y,n=length(y))
               y[.isEqual01(y)&y<0.4] <- 0
               Y <- scaleX.fct(y)
               yf <- prettyNum(y)
               Y <- scaleY.fct(y)
               i01 <- !.isEqual01(Y)
               yf <- yf[i01]
               Yi <- Y
               Y <- Y[i01]
               i0 <- any(!i01&Yi<0.5)
               i1 <- any(!i01&Yi>0.5)
               if(i0){ yf <- c(NA,yf); Y <- c(0, Y)}
               if(i1){ yf <- c(yf,NA); Y <- c(Y, 1)}
               axis(2,at=Y,labels=yf)
               if(i0) axis(2,at=0,labels=expression(-infinity))
               if(i1) axis(2,at=1,labels=expression(infinity))
            }else{
               if(is.null(ylim)){ ylim <- c(-Inf,Inf)}else{
                  if(is.na(ylim[1])) ylim[1] <- -Inf
                  if(is.na(ylim[2])) ylim[2] <- Inf }
               y.ticks <- sort(unique(y.ticks[!is.na(y.ticks)]))
               yf <- pmin(pmax(y.ticks[is.finite(y.ticks)],ylim[1]),ylim[2])
               Yf <- scaleY.fct(yf)
               axis(2,at=Yf,labels=yf)
               if(-Inf %in% y.ticks) axis(2,at=0,labels=expression(-infinity))
               if(Inf %in% y.ticks)  axis(2,at=1,labels=expression(infinity))
            }
            box()
        }else{
            if(!is.null(y.ticks)){
               if(is.null(ylim)){ ylim <- c(-Inf,Inf)}else{
                  if(is.na(ylim[1])) ylim[1] <- -Inf
                  if(is.na(ylim[2])) ylim[2] <- Inf }
               y.ticks <- sort(unique(y.ticks[!is.na(y.ticks)]))
               yf <- pmin(pmax(y.ticks[is.finite(y.ticks)],ylim[1]),ylim[2])
               axis(2,at=yf,labels=yf)
               if(-Inf %in% y.ticks) axis(2,at=0,labels=expression(-infinity))
               if(Inf %in% y.ticks)  axis(2,at=1,labels=expression(infinity))
               box()
           }
        }
   return(invisible(NULL))
}
if(FALSE){
  set.seed(1); x<- sort(c(-10,rnorm(100),10))
  xlim0 <- c(-2,1.6)
  ylim0 <- c(-0.8,1)
  xlim01 <- ex0(xlim0)
  ylim01 <- ex0(ylim0)
  xlim0 <- NULL
  ylim0 <- NULL
  xlim01 <- NULL
  ylim01 <-NULL
  distroptions("DistrResolution"=0.000001)
  res3 <- .rescalefct(x, fct=function(s) sin(s),
            scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T,
            xlim=xlim0,ylim=ylim0,dots=list(NULL))
  ex1 <- function(x)log(x/(1-x))
  ex0 <- function(x)exp(x)/(exp(x)+1)
  res4 <- .rescalefct(x, fct=function(s) sin(s),
            scaleX=T, scaleX.fct=ex0,
            scaleX.inv = ex1, scaleY=T,
            xlim=xlim0,ylim=ylim0,dots=list(NULL))
  plot(res3$X,res3$Y,axes=F, xlim=xlim01,ylim=ylim01)
  .plotRescaledAxis(scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T, scaleY.fct=pnorm,
            scaleY.inv = qnorm, X= res3$X, xlim=xlim0,ylim=ylim0, m = 19)
  plot(res3$X,res3$Y,axes=F, xlim=xlim01,ylim=ylim01)
  .plotRescaledAxis(scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T, scaleY.fct=pnorm,
            scaleY.inv = qnorm, X= res3$X, xlim=xlim0,ylim=ylim0,
            x.ticks = c(-100,-3,-1,-0.3,0,0.5,2,5,100),
            y.ticks = c(-1,-0.7,-0.1,0,0.2,.5,1))
  plot(res3$X,res3$Y,axes=F, xlim=xlim01,ylim=ylim01)
  .plotRescaledAxis(scaleX=T, scaleX.fct=function(x)exp(x)/(exp(x)+1),
            scaleX.inv = function(x)log(x/(1-x)), scaleY=T, scaleY.fct=pnorm,
            scaleY.inv = qnorm, X= res3$X, xlim=xlim0,ylim=ylim0,
            x.ticks = c(-Inf,-3,-1,-0.3,0,0.5,2,5,Inf),
            y.ticks = c(-1,-0.7,-0.1,0,0.2,.5,1))

}

.legendCoord <- function(x, scaleX, scaleX.fct, scaleY, scaleY.fct){
# rescaled legend coordinates axes acc. to logicals scaleX, scaleY
# return value: transformed legend coordinates
                if (is.character(x)) return(x)
                x1 <- if(scaleX) scaleX.fct(x[1]) else x[1]
                x2 <- if(scaleY) scaleY.fct(x[2]) else x[2]
                return(c(x1,x2))
            }
