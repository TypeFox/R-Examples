### R code from vignette source 'HistogramTools-quickref.Rnw'

###################################################
### code chunk number 1: HistogramTools-quickref.Rnw:28-32
###################################################
options(width=50)
library(HistogramTools)
set.seed(0)
ht.version <- packageDescription("HistogramTools")$Version


###################################################
### code chunk number 2: HistogramTools-quickref.Rnw:63-68 (eval = FALSE)
###################################################
## h <- hist(runif(100, 0, 100),
##           breaks=seq(from=0,to=200,by=5))
## plot(TrimHistogram(h))
## plot(SubsetHistogram(h, maxbreak=70))
## plot(MergeBuckets(h, adj.buckets=2))


###################################################
### code chunk number 3: HistogramTools-quickref.Rnw:71-73
###################################################
h <- hist(runif(100, 0, 100),
          breaks=seq(from=0,to=200,by=5), plot=F)


###################################################
### code chunk number 4: exhist
###################################################
# top, left, bottom, right
old.mar <- par(mar=c(3,4,1,1))
par(mfrow=c(2,2))
plot(h, main="Histogram h")
plot(TrimHistogram(h), main="TrimHistogram(h)")
plot(SubsetHistogram(h, max=70), main="SubsetHistogram(h, max=70)")
plot(MergeBuckets(h, 4), main="MergeBuckets(h, 4)")
par(mar=old.mar)


###################################################
### code chunk number 5: errorhist
###################################################
par(mfrow=c(1,2), par(mar=c(5,4,4,0)+0.1))
PlotEMDCC(h)
PlotKSDCC(h)
EMDCC(h)
KSDCC(h)


###################################################
### code chunk number 6: HistogramTools-quickref.Rnw:135-137 (eval = FALSE)
###################################################
## hist.msg <- as.Message(h)
## length(hist.msg$serialize(NULL))


###################################################
### code chunk number 7: HistogramTools-quickref.Rnw:139-143
###################################################
if(require(RProtoBuf)) {
  hist.msg <- as.Message(h)
  length(hist.msg$serialize(NULL))
}


