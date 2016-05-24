### R code from vignette source 'feature.Rnw'

###################################################
### code chunk number 1: feature.Rnw:67-73 (eval = FALSE)
###################################################
## library(feature)
## data(earthquake)
## eq3 <- log10(-earthquake[,3])
## eq3.fs <- featureSignif(eq3, bw=0.1)
## plot(eq3.fs, xlab="-log(-depth)")
## eq3.SiZer <- SiZer(eq3, bw=c(0.05, 0.5), xlab="-log(-depth)")


###################################################
### code chunk number 2: feature.Rnw:76-87
###################################################
library(feature)
data(earthquake)
eq3 <- log10(-earthquake[,3])
par(mar=c(4.5,4,2,2))
layout(matrix(1:2, nrow=2))
eq3.fs <- featureSignif(eq3, bw=0.1)
plot(eq3.fs, xlab="-log(-depth)", addSignifGradRegion=TRUE, addData=TRUE)
xlim <- par()$usr[1:2]
eq3.SiZer <- SiZer(eq3, xlim=xlim, bw=c(0.05, 0.5), xlab="-log(-depth)")
abline(h=log(0.1))
layout(1)


###################################################
### code chunk number 3: feature.Rnw:103-107 (eval = FALSE)
###################################################
## library(MASS)
## data(geyser)
## geyser.fs <- featureSignif(geyser, bw=c(4.5, 0.37))
## plot(geyser.fs, addSignifCurvRegion=TRUE)


###################################################
### code chunk number 4: feature.Rnw:109-113
###################################################
library(MASS)
data(geyser)
geyser.fs <- featureSignif(geyser, bw=c(4.5, 0.37), gridsize=c(51,51))
plot(geyser.fs, addSignifCurvRegion=TRUE)


###################################################
### code chunk number 5: feature.Rnw:124-125 (eval = FALSE)
###################################################
## plot(geyser.fs, addSignifCurvData=TRUE)


###################################################
### code chunk number 6: feature.Rnw:127-128
###################################################
plot(geyser.fs, addSignifCurvData=TRUE)


###################################################
### code chunk number 7: feature.Rnw:136-137
###################################################
names(geyser.fs)


