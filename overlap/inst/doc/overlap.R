### R code from vignette source 'overlap.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(continue="  ")


###################################################
### code chunk number 2: loadData
###################################################
library(overlap)
data(kerinci)
head(kerinci)
table(kerinci$Zone)
summary(kerinci$Sps)
range(kerinci$Time)


###################################################
### code chunk number 3: convertToRadians
###################################################
timeRad <- kerinci$Time * 2 * pi


###################################################
### code chunk number 4: singleDensityCurve
###################################################
tig2 <- timeRad[kerinci$Zone == 2 & kerinci$Sps == 'tiger']
densityPlot(tig2, rug=TRUE)


###################################################
### code chunk number 5: smoothing
###################################################
par(mfrow=2:1)
densityPlot(tig2, rug=TRUE, adjust=2)
text(2, 0.07, "adjust = 2")
densityPlot(tig2, rug=TRUE, adjust=0.2)
text(2, 0.1, "adjust = 0.2")


###################################################
### code chunk number 6: tigerMacaque2
###################################################
tig2 <- timeRad[kerinci$Zone == 2 & kerinci$Sps == 'tiger']
mac2 <- timeRad[kerinci$Zone == 2 & kerinci$Sps == 'macaque']
min(length(tig2), length(mac2))
tigmac2est <- overlapEst(tig2, mac2)
tigmac2est 
overlapPlot(tig2, mac2, main="Zone 2")
legend('topright', c("Tigers", "Macaques"), lty=c(1,2), col=c(1,4), bty='n')


###################################################
### code chunk number 7: bootstrap1
###################################################
tig2boot <- resample(tig2, 1000)
dim(tig2boot)
mac2boot <- resample(mac2, 1000)
dim(mac2boot)


###################################################
### code chunk number 8: bootstrap2
###################################################
tigmac2 <- bootEst(tig2boot, mac2boot, adjust = c(NA, 1, NA))  # takes a few seconds
dim(tigmac2)
BSmean <- colMeans(tigmac2)
BSmean


###################################################
### code chunk number 9: bootstrapCI
###################################################
tmp <- tigmac2[, 2] # Extract the required column of the matrix
bootCI(tigmac2est[2], tmp)


###################################################
### code chunk number 10: bootstrapCIlogit
###################################################
tmp <- tigmac2[, 2] # Extract the required column of the matrix
bootCIlogit(tigmac2est[2], tmp)


