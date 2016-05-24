library(diveMove)

zz <- system.file(file.path("data", "dives.csv"),
                  package="diveMove", mustWork=TRUE)
(sealX <- readTDR(zz, concurrentCols=4:6, speed=TRUE,
                  sep=";", na.strings="", as.is=TRUE))
(dcalib <- calibrateDepth(sealX, dry.thr=3610, zoc.method="offset", offset=3))
(dcalib <- calibrateDepth(sealX, zoc.method="offset", offset=3,
                          smooth.par=0.1, knot.factor=3))

###_+ Check all calibrateDepth() procedure

###_ : Check phase detection
detp <- diveMove:::.detPhase(getTime(sealX), getDepth(sealX), dry.thr=70,
                             wet.thr=3610, interval=getDtime(sealX))
###_ : Check zoc
zd <- diveMove:::.zoc(getTime(sealX), getDepth(sealX),
                      method="offset", control=list(offset=3))
if (!is.null(zd)) sealX@depth <- zd
###_ : Check dive detection
detd <- diveMove:::.detDive(getDepth(sealX), detp[[2]], 4)
###_ : Check labelling of dive phases
phaselabs <- diveMove:::.labDivePhase(sealX, detd[, 1], smooth.par=0.1,
                                      knot.factor=3, descent.crit.q=0,
                                      ascent.crit.q=0.1)

vcalib <- calibrateSpeed(dcalib, z=0, cex.pts=0.2)


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
