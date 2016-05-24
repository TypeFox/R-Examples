### R code from vignette source 'blockcluster_tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(blockcluster)
bc.version <- packageDescription("blockcluster")$Version
bc.date <- packageDescription("blockcluster")$Date


###################################################
### code chunk number 2: blockcluster_tutorial.Rnw:408-410
###################################################
defaultstrategy <- coclusterStrategy()
summary(defaultstrategy)


###################################################
### code chunk number 3: blockcluster_tutorial.Rnw:418-419
###################################################
newstrategy <- coclusterStrategy(nbtry=5, nbxem=10, algo='BCEM')


###################################################
### code chunk number 4: blockcluster_tutorial.Rnw:587-590
###################################################
data("binarydata")
out<-coclusterBinary(binarydata, nbcocluster=c(2,3))
summary(out)


###################################################
### code chunk number 5: blockcluster_tutorial.Rnw:600-601
###################################################
plot(out, asp = 0)


###################################################
### code chunk number 6: blockcluster_tutorial.Rnw:607-608
###################################################
plot(out, type = 'distribution')


