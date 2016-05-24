### R code from vignette source 'BigData.rnw'

###################################################
### code chunk number 1: BigData.rnw:3-4
###################################################
options(width=50)


###################################################
### code chunk number 2: BigData.rnw:56-61
###################################################
library(RecordLinkage)

showClass("RLBigData")
showClass("RLBigDataDedup")
showClass("RLBigDataLinkage")


###################################################
### code chunk number 3: BigData.rnw:71-84
###################################################
# deduplicate dataset with two blocking iterations and string comparison
data(RLdata500)
data(RLdata10000)
rpairs1 <- RLBigDataDedup(RLdata500, identity = identity.RLdata500, blockfld = list(1,3),
  strcmp = 1:4)

# link two datasets with phonetic code, exclude lname_c2
s1 <- 471:500
s2 <- sample(1:10000, 300)
identity2 <- c(identity.RLdata500[s1], rep(NaN, length(s2)))
dataset <- rbind(RLdata500[s1,], RLdata10000[s2,])
rpairs2 <- RLBigDataLinkage(RLdata500, dataset, identity1 = identity.RLdata500,
  identity2 = identity2, phonetic = 1:4, exclude = "lname_c2")


###################################################
### code chunk number 4: BigData.rnw:95-100
###################################################
train <- getMinimalTrain(compare.dedup(RLdata500, identity = identity.RLdata500,
  blockfld = list(1,3)))
rpairs1 <- RLBigDataDedup(RLdata500, identity = identity.RLdata500)
classif <- trainSupv(train, "rpart", minsplit=2)
result <- classifySupv(classif, rpairs1)


###################################################
### code chunk number 5: BigData.rnw:105-106
###################################################
showClass("RLResult")


###################################################
### code chunk number 6: BigData.rnw:110-112
###################################################
getTable(result)
getErrorMeasures(result)


###################################################
### code chunk number 7: BigData.rnw:123-126
###################################################
rpairs1 <- epiWeights(rpairs1)
result <- epiClassify(rpairs1, 0.5)
getTable(result)


###################################################
### code chunk number 8: BigData.rnw:140-141
###################################################
getPairs(result, min.weight=0.7, filter.link="link")


###################################################
### code chunk number 9: BigData.rnw:148-151
###################################################
getFalsePos(result)
getFalseNeg(result)



