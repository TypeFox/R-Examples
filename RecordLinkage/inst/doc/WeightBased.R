### R code from vignette source 'WeightBased.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: WeightBased.rnw:3-4
###################################################
options(width=50)


###################################################
### code chunk number 2: WeightBased.rnw:25-26
###################################################
library(RecordLinkage)


###################################################
### code chunk number 3: WeightBased.rnw:33-35
###################################################
data(RLdata500)
RLdata500[1:5,]


###################################################
### code chunk number 4: WeightBased.rnw:44-46
###################################################
pairs=compare.dedup(RLdata500,identity=identity.RLdata500,blockfld=list(c(5,6),c(6,7),c(5,7)))
summary(pairs)


###################################################
### code chunk number 5: WeightBased.rnw:54-56
###################################################
pairs=emWeights(pairs)
hist(pairs$Wdata, plot=FALSE)


###################################################
### code chunk number 6: WeightBased.rnw:69-70
###################################################
getPairs(pairs,30,20)


###################################################
### code chunk number 7: WeightBased.rnw:72-73
###################################################
getPairs(pairs,30,20)[23:36,]


###################################################
### code chunk number 8: WeightBased.rnw:75-77
###################################################
pairs=emClassify(pairs, threshold.upper=24, threshold.lower=-7)
summary(pairs)


###################################################
### code chunk number 9: WeightBased.rnw:85-91
###################################################
possibles <- getPairs(pairs, show="possible")
possibles[1:6,]
links=getPairs(pairs,show="links", single.rows=TRUE)
link_ids <- links[, c("id1", "id2")]
link_ids



