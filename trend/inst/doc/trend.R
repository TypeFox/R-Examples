### R code from vignette source 'trend.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: trend.Rnw:85-89
###################################################
require(trend)
data(maxau)
Q <- maxau[,"Q"]
mk.test(Q)


###################################################
### code chunk number 2: trend.Rnw:109-111
###################################################
require(trend)
smk.test(nottem)


###################################################
### code chunk number 3: trend.Rnw:119-121
###################################################
require(trend)
csmk.test(nottem)


###################################################
### code chunk number 4: trend.Rnw:131-134
###################################################
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
cor.test(s,Q, meth="spearman")


###################################################
### code chunk number 5: trend.Rnw:139-143
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.mk.test(s,Q)


###################################################
### code chunk number 6: trend.Rnw:152-156
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.cor.trend.test(s,Q, "spearman")


###################################################
### code chunk number 7: trend.Rnw:187-190
###################################################
require(trend)
s <- maxau[,"s"]
sens.slope(s)


###################################################
### code chunk number 8: trend.Rnw:205-207
###################################################
require(trend)
sea.sens.slope(nottem)


###################################################
### code chunk number 9: trend.Rnw:234-237
###################################################
require(trend)
data(PagesData)
pettitt.test(PagesData)


