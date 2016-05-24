### R code from vignette source 'Supervised.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Supervised.rnw:3-4
###################################################
options(width=50)


###################################################
### code chunk number 2: Supervised.rnw:27-28
###################################################
library(RecordLinkage)


###################################################
### code chunk number 3: Supervised.rnw:35-41
###################################################
data(RLdata500)
data(RLdata10000)
train_pairs=compare.dedup(RLdata10000, identity=identity.RLdata10000,
  n_match=500, n_non_match=500)

eval_pairs=compare.dedup(RLdata500,identity=identity.RLdata500)


###################################################
### code chunk number 4: Supervised.rnw:52-55
###################################################
model_rpart=trainSupv(train_pairs, method="rpart")
model_bagging=trainSupv(train_pairs, method="bagging")
model_svm=trainSupv(train_pairs, method="svm")


###################################################
### code chunk number 5: Supervised.rnw:64-67
###################################################
result_rpart=classifySupv(model_rpart, eval_pairs)
result_bagging=classifySupv(model_bagging, eval_pairs)
result_svm=classifySupv(model_svm, eval_pairs)


###################################################
### code chunk number 6: Supervised.rnw:73-74
###################################################
texSummary(result_rpart)


###################################################
### code chunk number 7: Supervised.rnw:78-79
###################################################
texSummary(result_bagging)


###################################################
### code chunk number 8: Supervised.rnw:83-84
###################################################
texSummary(result_svm)


