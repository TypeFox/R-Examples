## ----setup, include=FALSE------------------------------------------------
# fig.width=5.5, 
require(knitr)
opts_chunk$set(fig.align='center', fig.height=5, warning=FALSE,
               dev='pdf', prompt=TRUE, comment=NA, highlight=FALSE, tidy=FALSE)

## ----nplr, message=FALSE, warning=FALSE----------------------------------
require(nplr)

## ----test1---------------------------------------------------------------
path <- system.file("extdata", "pc3.txt", package="nplr")
pc3 <- read.delim(path)
np1 <- nplr(x=pc3$CONC, y=pc3$GIPROP)

## ----out1----------------------------------------------------------------
np1

## ----plot1, include=FALSE------------------------------------------------
plot(np1, cex.main = 1.2,
     main="PC-3 cell line. Response to Thioguanine")

## ----fig1, ref.label='plot1'---------------------------------------------
plot(np1, cex.main = 1.2,
     main="PC-3 cell line. Response to Thioguanine")

## ----custom, include=FALSE-----------------------------------------------
op <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=.5, showInfl=TRUE,
     main="Default 'nplr' plot", cex.main=1.5)
x1 <- getX(np1); y1 <- getY(np1)
x2 <- getXcurve(np1); y2 <- getYcurve(np1)
plot(x1, y1, pch=15, cex=2, col="tan1", xlab=expression(Log[10](conc)),
     ylab="Prop", main="Custom plot", cex.main=1.5)
lines(x2, y2, lwd=5, col="seagreen4")
par(op)

## ----plotCustom, ref.label='custom', fig.width=12, fig.height=5.5, echo=3:7----
op <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(np1, pcol="grey40", lcol="skyblue1", showEstim=.5, showInfl=TRUE,
     main="Default 'nplr' plot", cex.main=1.5)
x1 <- getX(np1); y1 <- getY(np1)
x2 <- getXcurve(np1); y2 <- getYcurve(np1)
plot(x1, y1, pch=15, cex=2, col="tan1", xlab=expression(Log[10](conc)),
     ylab="Prop", main="Custom plot", cex.main=1.5)
lines(x2, y2, lwd=5, col="seagreen4")
par(op)

## ----getGoodness1--------------------------------------------------------
getGoodness(np1)

## ----getStdErr1----------------------------------------------------------
getStdErr(np1)

## ----getPar1-------------------------------------------------------------
getPar(np1)

## ----getAUC1-------------------------------------------------------------
getAUC(np1)

## ----getEstimates1-------------------------------------------------------
getEstimates(np1)

## ----getTarget1----------------------------------------------------------
getEstimates(np1, .5)

## ----getTargetCI---------------------------------------------------------
getEstimates(np1, c(.25, .5, .75), conf.level=.90)

## ----test2---------------------------------------------------------------
path <- system.file("extdata", "mcf7.txt", package="nplr")
mcf7 <- read.delim(path)
np2 <- nplr(x=mcf7$CONC, y=mcf7$GIPROP)

## ----plot2---------------------------------------------------------------
plot(np2, showSDerr = TRUE, lwd = 4 , cex.main=1.25,
     main="Cell line MCF-7. Response to Irinotecan")

## ----testWeights, message=FALSE, warning=FALSE---------------------------
x <- mcf7$CONC
y <- mcf7$GIPROP
noweight <- nplr(x, y, LPweight=0, silent=TRUE)
sdw <- nplr(x, y, method="sdw", silent=TRUE)
gw <-  nplr(x, y, method="sdw", LPweight=1.5, silent=TRUE)

## ----plotWeights, fig.width=12, fig.height=10, echo=2:5------------------
par(mfrow=c(2,2))
plot(np2, showEstim=.5, main="residuals weights")
plot(noweight, showEstim=.5, main="No weight")
plot(sdw, showEstim=.5, main="Stdev weights")
plot(noweight, showEstim=.5, main="general weights")
par(op)

## ----loadProg------------------------------------------------------------
path <- system.file("extdata", "prog.txt", package="nplr")
prog <- read.delim(path)

## ----test3---------------------------------------------------------------
x <- prog$time
yp <- convertToProp(prog$prog, T0 = 5, Ctrl = 102)
np3 <- nplr(x, yp, useLog=FALSE)

## ----getInf3-------------------------------------------------------------
getInflexion(np3)

## ----plot3---------------------------------------------------------------
plot(np3, showInfl=TRUE, xlab="Time (hrs)", cex.main=1.5, cex.lab=1.2,
     ylab="Prop. of control", main="Progression")

## ----npar----------------------------------------------------------------
plot(x, yp , cex.main=1.5, cex.lab=1.2,
     main="The n-parameter effect", xlab="Time", ylab="Progression")
le <- c()
for(i in 2:5){
  test <- nplr(x, yp, npars=i, useLog=FALSE)
  lines(getXcurve(test), getYcurve(test), lwd=2, col=i)
  points(getInflexion(test), pch=19, cex=1.25, col=i)
  gof <- getGoodness(test)
  le <- c(le, sprintf("%s-P: GOF=%s", i, round(gof, 4)))
}
legend("bottomright", legend=le, lwd=2, pch=19, col=2:5, bty="n")

## ----overlay-------------------------------------------------------------
path <- system.file("extdata", "multicell.tsv", package="nplr")
multicell <- read.delim(path)

# Computing models (to store in a list)
cellsList <- split(multicell, multicell$cell)
Models <- lapply(cellsList, function(tmp){
  nplr(tmp$conc, tmp$resp, silent = TRUE)
  })

# Visualizing
overlay(Models, xlab = expression(Log[10](Conc.)), ylab = "Resp.",
  main="Superimposing multiple curves", cex.main=1.5, lwd = 3)

