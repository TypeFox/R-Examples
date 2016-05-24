### R code from vignette source 'diveMove.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: R-opts
###################################################
options(width=34, digits=4)


###################################################
### code chunk number 2: startup
###################################################
library(diveMove)


###################################################
### code chunk number 3: diveMove.Rnw:164-165 (eval = FALSE)
###################################################
## example(diveMove)


###################################################
### code chunk number 4: dives-con
###################################################
fp <- file.path("data", "dives.csv")
sfp <- system.file(fp, package="diveMove")


###################################################
### code chunk number 5: readin-csv
###################################################
srcfn <- basename(sfp)
tdrXcsv <- read.csv(sfp, sep=";")


###################################################
### code chunk number 6: create-tdr
###################################################
ddtt.str <- paste(tdrXcsv$date, tdrXcsv$time)
ddtt <- strptime(ddtt.str, 
                 format="%d/%m/%Y %H:%M:%S")
time.posixct <- as.POSIXct(ddtt, tz="GMT")
tdrX <- createTDR(time=time.posixct, 
                  depth=tdrXcsv$depth,
                  concurrentData=tdrXcsv[, -c(1:3)], 
                  dtime=5, file=srcfn)
## Or a TDRspeed object, since we know we have
## speed measurements:
tdrX <- createTDR(time=time.posixct, 
                  depth=tdrXcsv$depth,
                  concurrentData=tdrXcsv[, -c(1:3)], 
                  dtime=5, file=srcfn, 
                  speed=TRUE)


###################################################
### code chunk number 7: readin-tdr (eval = FALSE)
###################################################
## fp <- file.path("data", "dives.csv")
## sfp <- system.file(fp, package="diveMove")
## tdrX <- readTDR(sfp, speed=TRUE, sep=";", 
##                 na.strings="", as.is=TRUE)
## plotTDR(tdrX)


###################################################
### code chunk number 8: plot-tdr
###################################################
fp <- file.path("data", "dives.csv")
sfp <- system.file(fp, package="diveMove")
tdrX <- readTDR(sfp, speed=TRUE, sep=";", 
                na.strings="", as.is=TRUE)
plotTDR(tdrX, interact=FALSE, cex.lab=1.3)


###################################################
### code chunk number 9: diveMove.Rnw:333-334 (eval = FALSE)
###################################################
## dcalib <- calibrateDepth(tdrX)


###################################################
### code chunk number 10: diveMove.Rnw:356-359 (eval = FALSE)
###################################################
## dcalib <- calibrateDepth(tdrX, 
##                          zoc.method="offset",
##                          offset=3)


###################################################
### code chunk number 11: diveMove.Rnw:381-386 (eval = FALSE)
###################################################
## dcalib <- calibrateDepth(tdrX, 
##                          zoc.method="filter",
##                          k=c(3, 5760), 
##                          probs=c(0.5, 0.02), 
##                          na.rm=TRUE)


###################################################
### code chunk number 12: zoc
###################################################
dcalib <- calibrateDepth(tdrX, dive.thr=3, 
                         zoc.method="offset",
                         offset=3, descent.crit.q=0.01, 
                         ascent.crit.q=0,
                         knot.factor=20)


###################################################
### code chunk number 13: plot-gross-activity (eval = FALSE)
###################################################
## plotTDR(dcalib, concurVars=c("speed", "light"),
##         surface=TRUE)


###################################################
### code chunk number 14: plot-tdrcalibrate
###################################################
plotTDR(dcalib, concurVars=c("speed", "light"),
        surface=TRUE, interact=FALSE, cex.lab=1.3)


###################################################
### code chunk number 15: plot-dive-activity (eval = FALSE)
###################################################
## plotTDR(dcalib, diveNo=2:8, what="phases")


###################################################
### code chunk number 16: plot-tdr-dives
###################################################
plotTDR(dcalib, diveNo=2:8, what="phases", 
        interact=FALSE, cex.lab=1.3, 
        depth.lim=c(0, 80))


###################################################
### code chunk number 17: extract-dive (eval = FALSE)
###################################################
## extractDive(dcalib, diveNo=2:8)


###################################################
### code chunk number 18: tdr-extract
###################################################
getTDR(dcalib)


###################################################
### code chunk number 19: grossact1 (eval = FALSE)
###################################################
## getGAct(dcalib)


###################################################
### code chunk number 20: grossact2 (eval = FALSE)
###################################################
## getGAct(dcalib, "phase.id")


###################################################
### code chunk number 21: diveact-1 (eval = FALSE)
###################################################
## getDAct(dcalib)


###################################################
### code chunk number 22: dphaselab1 (eval = FALSE)
###################################################
## getDPhaseLab(dcalib)
## getDPhaseLab(dcalib, 20)


###################################################
### code chunk number 23: dphaselab2
###################################################
dphases <- getDPhaseLab(dcalib, c(100:300))


###################################################
### code chunk number 24: diveModel (eval = FALSE)
###################################################
## plotDiveModel(dcalib, diveNo=260)


###################################################
### code chunk number 25: plot-dive-model
###################################################
plotDiveModel(dcalib, diveNo=260)


###################################################
### code chunk number 26: extractdive
###################################################
sealX <- extractDive(dcalib, diveNo=c(100:300))
sealX


###################################################
### code chunk number 27: plot-phases (eval = FALSE)
###################################################
## plotTDR(sealX, phaseCol=dphases)


###################################################
### code chunk number 28: diveMove.Rnw:670-671
###################################################
options(width=105)


###################################################
### code chunk number 29: dive-summaries
###################################################
tdrXSumm1 <- diveStats(dcalib)
names(tdrXSumm1)
tbudget <- timeBudget(dcalib, ignoreZ=TRUE)
head(tbudget, 4)
trip.labs <- stampDive(dcalib, ignoreZ=TRUE)
tdrXSumm2 <- data.frame(trip.labs, tdrXSumm1)
names(tdrXSumm2)


###################################################
### code chunk number 30: diveMove.Rnw:680-681
###################################################
options(width=34, digits=4)


###################################################
### code chunk number 31: calibrate-speed (eval = FALSE)
###################################################
## vcalib <- calibrateSpeed(dcalib, tau=0.1, 
##                          contour.level=0.1,
##                          z=1, bad=c(0, 0), 
##                          cex.pts=0.2)


###################################################
### code chunk number 32: plot-speed-calibration
###################################################
calibrateSpeed(dcalib, tau=0.1, 
               contour.level=0.1,
               z=1, bad=c(0, 0), 
               cex.pts=0.2)


