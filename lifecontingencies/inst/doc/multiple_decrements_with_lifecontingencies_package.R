### R code from vignette source 'multiple_decrements_with_lifecontingencies_package.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: multiple_decrements_with_lifecontingencies_package.Rnw:53-54
###################################################
options(width=80, prompt='R> ')


###################################################
### code chunk number 2: load
###################################################
library(lifecontingencies)


###################################################
### code chunk number 3: mdt1
###################################################
valdezDf<-data.frame(
		x=c(50:54),
		lx=c(4832555,4821937,4810206,4797185,4782737),
		hearth=c(5168, 5363, 5618, 5929, 6277),
		accidents=c(1157, 1206, 1443, 1679,2152),
		other=c(4293,5162,5960,6840,7631)
)
valdezMdt<-new("mdt",name="ValdezExample",table=valdezDf)


###################################################
### code chunk number 4: md3a (eval = FALSE)
###################################################
## print(valdezMdt)


###################################################
### code chunk number 5: md3b
###################################################
valdezDf<-as(valdezMdt,"data.frame")
require(markovchain)
valdezMarkovChainList<-as(valdezMdt,"markovchainList")


###################################################
### code chunk number 6: mdt4
###################################################
getOmega(valdezMdt)
getDecrements(valdezMdt)


###################################################
### code chunk number 7: summary
###################################################
summary(valdezMdt)


###################################################
### code chunk number 8: dx1
###################################################
dxt(valdezMdt,x=51,decrement="other")
dxt(valdezMdt,x=51,t=2, decrement="other")
dxt(valdezMdt,x=51)


###################################################
### code chunk number 9: dx2
###################################################
dxt(valdezMdt,x=51,t=2, decrement="other")
pxt(valdezMdt,x=50,t=3)
qxt(valdezMdt,x=53,t=2,decrement=1)


###################################################
### code chunk number 10: randomSamples
###################################################
rmdt(n = 2,object = valdezMdt,x = 50,t = 2)


###################################################
### code chunk number 11: act1
###################################################
myTable<-data.frame(x=c(16,17,18),
  lx=c(20000,17600,14520),
  da=c(1300,1870,2380),
  doc=c(1100,1210,1331)
)
myMdt<-new("mdt",table=myTable,name="Sample")


###################################################
### code chunk number 12: act2
###################################################
Axn.mdt(object=myMdt,x=16,i=.1,decrement="da")


