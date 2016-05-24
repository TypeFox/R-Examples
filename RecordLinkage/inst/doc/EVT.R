### R code from vignette source 'EVT.rnw'

###################################################
### code chunk number 1: EVT.rnw:3-4
###################################################
options(width=50)


###################################################
### code chunk number 2: EVT.rnw:28-29
###################################################
library(RecordLinkage)


###################################################
### code chunk number 3: EVT.rnw:32-37
###################################################
data(RLdata500)
rpairs=compare.dedup(RLdata500,identity=identity.RLdata500,
  blockfld=list(1,3,5,6,7),strcmp=1:4)
rpairs=emWeights(rpairs)



###################################################
### code chunk number 4: EVT.rnw:54-55
###################################################
## Not run: getParetoThreshold(rpairs)


###################################################
### code chunk number 5: EVT.rnw:59-60
###################################################
plotMRL(rpairs)


###################################################
### code chunk number 6: EVT.rnw:68-73
###################################################
plotMRL(rpairs)
abline(v=c(1.2,12.8),col="red",lty="dashed")
l=mrl(rpairs$Wdata)
range=l$x>1.2 & l$x < 12.8
points(l$x[range], l$y[range],col="red",type="l")


###################################################
### code chunk number 7: EVT.rnw:83-86
###################################################
threshold=getParetoThreshold(rpairs,interval=c(1.2,12.8))
result=emClassify(rpairs,threshold)
summary(result)


