## ---- echo=TRUE, message=FALSE-------------------------------------------
library(secr)

## ---- echo=FALSE, eval=TRUE--------------------------------------------------------
options(digits=6, width=85)

## ---- echo=FALSE, eval=FALSE-------------------------------------------------------
#  par(mfrow=c(1,1), pty='m', mar=c(4,6,2,6), las=1, bty='l',
#      xpd=T, cex=1.2, xaxs='i', yaxs='i')
#  plot(0,0, type='n', xlim=c(0,800), ylim=c(0,0.15),
#     xlab='', ylab='', lwd=2, col='blue', axes = FALSE)
#  plot(stoat.model.HN, limits=FALSE, xv=0:800, ylim=c(0,0.12),
#     xlab='', ylab='', lwd=2, col='blue', add = TRUE)
#  axis (1)
#  mtext(side=1, line=2.5, 'Distance  (m)', cex=1.2)
#  axis (2, at=c(0, 0.05, 0.10, 0.15))
#  mtext(side=2, line=3.5, 'Detection probability', cex=1.2, las=0)
#  plot(stoat.model.EX, add=T, limits=F, xv=1:800, col='green', lwd=2)
#  legend (500,0.12, lwd=2, col=c('blue', 'green'),
#      legend=c('halfnormal','exponential'), bty='n')

## ---- eval=FALSE-------------------------------------------------------------------
#  RShowDoc ("secr-manual", package = "secr")

## ---- eval=FALSE-------------------------------------------------------------------
#  news (package = "secr")

## ---- eval=FALSE-------------------------------------------------------------------
#  secr.fit(captdata, model = g0~t)

## ---- echo=FALSE, eval=TRUE--------------------------------------------------------
options(digits = 6, width = 85)       
library(secr)                                               # load package                             
myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')    # import data using XY format

## ---- eval=FALSE-------------------------------------------------------------------
#  library(secr)                                               # load package
#  oldwd <- setwd(system.file('extdata', package = 'secr'))    # change working folder
#  myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')    # import data using XY format
#  setwd(oldwd)                                                # reset working folder

