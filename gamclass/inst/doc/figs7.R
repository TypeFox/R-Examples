## ----setup, cache=FALSE, echo=FALSE-----------------------------------
library(knitr)
options(replace.assign=FALSE,width=72)
opts_chunk$set(fig.path='figs/gam-', cache.path='cache/gam-',
               fig.align='center', dev='pdf', fig.width=3.5,
               fig.height=3.5, fig.show='hold', pars=TRUE,
               tidy=FALSE,  comment=NA)
knit_hooks$set(pars=function(before, options, envir){
if (before && options$fig.show!='none') par(mar=c(4,4,1.6,.1),
              font.main=1, cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),
              tcl=-.3)
              par(options$pars)
}, crop=hook_pdfcrop)
pdf.options(pointsize=12)
oldopt <- options(digits=4)

## ----fig7_1, eval=TRUE, echo=TRUE-------------------------------------
fig7.1 <- function(){
Erie <- greatLakes[,"Erie"]
plot(Erie, xlab="",
     ylab="Level (m)")
}

## ----fig7_2, eval=TRUE, echo=TRUE-------------------------------------
fig7.2 <- function(){
    Erie <- greatLakes[,"Erie"]
    opar <- par(oma=c(0,0,4,0))
    lag.plot(Erie, lags=3,
             do.lines=FALSE,
             layout=c(2,3), main="")
    mtext(side=3, line=3, adj=-0.155,
          "A: Lag plots, for lags 1, 2 and 3 respectively", cex=1)
    par(fig=c(0,1,0,0.6), new=TRUE)
    par(mar=c(2.75, 3.1, 3.6, 1.6))
    acf(Erie, main="", xlab="")
    mtext(side=3, line=0.5, "B: Autocorrelation estimates at successive lags",
          adj=-0.35, cex=1)
    mtext(side=1, line=1.75, "Lag", cex=1)
    par(fig=c(0,1,0,1))
    par(opar)
}

## ----fig7_3, eval=TRUE, echo=TRUE-------------------------------------
fig7.3 <- function(){
    Erie <- greatLakes[,"Erie"]
    df <-  data.frame(height=as.vector(Erie), year=time(Erie))
    obj <- gam(height ~ s(year), data=df)
    plot(obj, shift=mean(df$height), residuals=T, pch=1, xlab="")
}

## ----fig7_4, eval=TRUE, echo=TRUE-------------------------------------
fig7.4 <- function(){
    if(!require(forecast))return("Package 'forecast' must be installed")
    Erie <- greatLakes[,"Erie"]
    assign('Erie', Erie, pos=1)
    erie.ar <- ar(Erie)
    plot(forecast(erie.ar, h=15), ylab="Lake level (m)")
}

## ----fig7_5, eval=TRUE, echo=TRUE-------------------------------------
fig7.5 <- function(mf=3,nf=2){
    opar <- par(mfrow=c(mf,nf), mar=c(0.25, 4.1, 0.25, 1.1))
    npanel <- mf*nf
    for(i in 1:npanel){
        df <- data.frame(x=1:200, y=arima.sim(list(ar=0.7), n=200))
        df.gam <- gam(y ~ s(x), data=df)
        plot(df.gam, residuals=TRUE)
    }
    par(opar)
}

## ----fig7_6, eval=TRUE, echo=TRUE-------------------------------------
fig7.6 <- function(){
    mdbRain.gam <- gam(mdbRain ~ s(Year) + s(SOI), data=bomregions2012)
    plot(mdbRain.gam, residuals=TRUE, se=2, pch=1, cex=0.5, select=1)
    plot(mdbRain.gam, residuals=TRUE, se=2, pch=1, cex=0.5, select=2)
}

## ----figs7-do, eval=TRUE, message=FALSE, warning=FALSE----------------
pkgs <- c("DAAG","mgcv","splines","forecast")
z <- sapply(pkgs, require, character.only=TRUE, warn.conflicts=FALSE)
if(any(!z)){
  notAvail <- paste(names(z)[!z], collapse=", ")
  print(paste("The following packages require to be installed:", notAvail))
}

## ----fig7_1x, eval=TRUE, echo=TRUE, fig.width=4.5, fig.height=2.5, out.width="0.75\\textwidth"----
fig7.1()

## ----fig7_2x, eval=TRUE, echo=TRUE, fig.width=6, fig.height=4.5, out.width="0.75\\textwidth"----
fig7.2()

## ----fig7_3x, eval=TRUE, echo=TRUE, fig.width=6, fig.height=3.5, out.width="0.75\\textwidth"----
fig7.3()

## ----fig7_4x, eval=TRUE, echo=TRUE, fig.width=4.5, fig.height=2.25, out.width="0.75\\textwidth"----
fig7.4()

## ----fig7_5x, eval=TRUE, echo=TRUE, fig.width=6, fig.height=3.75, out.width="0.75\\textwidth"----
fig7.5()

## ----fig7_6x, eval=TRUE, echo=TRUE, fig.width=6, fig.height=3.25, pars=list(mfrow=c(1,2)), out.width="0.8\\textwidth"----
fig7.6()

