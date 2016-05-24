## ----eval = FALSE--------------------------------------------------------
#  read.capthist(captfile, trapfile, detector = "multi", fmt = c("trapID", "XY"),
#      noccasions = NULL, covnames = NULL, trapcovnames = NULL, cutval = NULL,
#      verify = TRUE, noncapt = "NONE", ...)

## ---- message=FALSE------------------------------------------------------
library(secr)
setwd (system.file('extdata', package = 'secr'))
stoatCH <- read.capthist('stoatcapt.txt', 'stoattrap.txt',
    detector = 'proximity')
summary(stoatCH)

## ---- eval=TRUE----------------------------------------------------------
write.capthist(signalCH, 'temp')  ## export data for demo
tempCH <- read.capthist('tempcapt.txt', 'temptrap.txt', detector = 'signal', cutval = 52.5)

## ---- eval = FALSE-------------------------------------------------------
#  read.capthist("captXY.txt", "perimeter.txt", fmt = 'XY', detector = "polygon")

## ---- eval = FALSE-------------------------------------------------------
#  temppoly <- read.traps(file = 'clipboard', detector = 'polygon')
#  tempcapt <- sim.capthist(temppoly, popn = list(D = 1, buffer = 1000), detectpar =
#                             list(g0 = 0.5, sigma = 250))
#  plot(tempcapt, label = TRUE, tracks = TRUE, title = 'Simulated detections within polygons')

## ------------------------------------------------------------------------
summary(subset(stoatCH, traps = 1:47, occasions = 1:5))

