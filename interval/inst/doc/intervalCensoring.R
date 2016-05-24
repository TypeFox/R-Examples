### R code from vignette source 'intervalCensoring.Rnw'

###################################################
### code chunk number 1: intervalCensoring.Rnw:102-103
###################################################
options(prompt = "R> ", continue = "+ ")


###################################################
### code chunk number 2: intervalCensoring.Rnw:716-720
###################################################
## set output line size
options(width = 65)
library(interval)
library(coin)


###################################################
### code chunk number 3: intervalCensoring.Rnw:723-726
###################################################
library("interval")
data("bcos", package = "interval")
head(bcos)


###################################################
### code chunk number 4: intervalCensoring.Rnw:736-738
###################################################
fit1<-icfit(Surv(left, right, type = "interval2")~treatment, data = bcos)
summary(fit1)


###################################################
### code chunk number 5: intervalCensoring.Rnw:759-760
###################################################
plot(fit1)


###################################################
### code chunk number 6: intervalCensoring.Rnw:775-777
###################################################
icout<-ictest(Surv(left, right, type = "interval2")~treatment, data = bcos)
icout


###################################################
### code chunk number 7: intervalCensoring.Rnw:788-789
###################################################
plot(fit1, dtype = "link")


###################################################
### code chunk number 8: intervalCensoring.Rnw:803-805
###################################################
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, 
initfit = icout$fit, scores = "logrank2")


###################################################
### code chunk number 9: intervalCensoring.Rnw:810-814
###################################################
L<-with(bcos, left)
R<-with(bcos, right)
trt<-with(bcos, treatment) 
ictest(L, R, trt, scores = "wmw", initfit = icout$fit)


###################################################
### code chunk number 10: intervalCensoring.Rnw:823-824
###################################################
set.seed(1232)


###################################################
### code chunk number 11: intervalCensoring.Rnw:826-828
###################################################
fakeTrtGrps<-sample(letters[1:4], nrow(bcos), replace = TRUE)
ictest(L, R, fakeTrtGrps)


###################################################
### code chunk number 12: intervalCensoring.Rnw:835-836
###################################################
set.seed(931)


###################################################
### code chunk number 13: intervalCensoring.Rnw:838-840
###################################################
fakeZ<-rnorm(nrow(bcos))
ictest(L, R, fakeZ, alternative = "less")


###################################################
### code chunk number 14: intervalCensoring.Rnw:847-848
###################################################
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, exact = TRUE, scores = "logrank1")


###################################################
### code chunk number 15: intervalCensoring.Rnw:863-864
###################################################
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, initfit = icout$fit, method = "scoretest", scores = "logrank2")


###################################################
### code chunk number 16: intervalCensoring.Rnw:871-873
###################################################
icoutHLY<-ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, initfit = icout$fit, method = "wsr.HLY",    mcontrol = mControl(nwsr = 99), scores = "logrank1")
icoutHLY


###################################################
### code chunk number 17: intervalCensoring.Rnw:883-884
###################################################
icoutHLY$fit$anypzero


###################################################
### code chunk number 18: intervalCensoring.Rnw:963-965
###################################################
data("ChickWeight", package = "datasets")
head(ChickWeight)


###################################################
### code chunk number 19: intervalCensoring.Rnw:977-978
###################################################
permTS(weight~Diet, data = ChickWeight, subset = Diet %in% c(3, 4) & Time == 21)


###################################################
### code chunk number 20: intervalCensoring.Rnw:981-984
###################################################
y3<-with(subset(ChickWeight, Time == 21 & Diet == 3), weight)
y4<-with(subset(ChickWeight, Time == 21 & Diet == 4), weight)
permTS(y3, y4)


###################################################
### code chunk number 21: intervalCensoring.Rnw:992-993
###################################################
permTS(y3[1:5], y4[1:5])


###################################################
### code chunk number 22: intervalCensoring.Rnw:1018-1019
###################################################
permKS(weight~Diet, data = ChickWeight, subset = Time == 21)


###################################################
### code chunk number 23: intervalCensoring.Rnw:1022-1025
###################################################
y<-ChickWeight[ChickWeight$Time == 21, "weight"]
g<-ChickWeight[ChickWeight$Time == 21, "Diet"]
permKS(y, g)


###################################################
### code chunk number 24: intervalCensoring.Rnw:1034-1035
###################################################
permTREND(y, as.numeric(g))


###################################################
### code chunk number 25: intervalCensoring.Rnw:1051-1053
###################################################
## set output line size
options(width = 65)


###################################################
### code chunk number 26: intervalCensoring.Rnw:1055-1059
###################################################
system.time(cm19c10<-chooseMatrix(length(y3)+length(y4),   length(y3)))
system.time(PC<-permControl(cm = cm19c10))
system.time(permTS(y3, y4, method = "exact.ce", control = PC))
system.time(permTS(y3, y4, method = "exact.network"))


###################################################
### code chunk number 27: intervalCensoring.Rnw:1176-1180
###################################################
icout<-ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, scores = "wmw")
wmw.scores<-icout$scores
logistic.scores<-ictest(Surv(left, right, type = "interval2")~treatment, data = bcos,    icFIT = icout$fit, scores = "general", dqfunc = function(x){ dlogis(qlogis(x))})$scores
max(abs(wmw.scores-logistic.scores))


###################################################
### code chunk number 28: intervalCensoring.Rnw:1213-1215
###################################################
library("coin")
independence_test(Surv(left, right, type = "interval2")~treatment, data = bcos, ytrafo = wlr_trafo)


###################################################
### code chunk number 29: intervalCensoring.Rnw:1220-1223
###################################################
SUBSET<-c(1:5, 50:65)
independence_test(Surv(left, right, type = "interval2")~treatment, data = bcos, subset = SUBSET, ytrafo = wlr_trafo, distribution = exact())
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, subset = SUBSET, method = "exact.network")


###################################################
### code chunk number 30: intervalCensoring.Rnw:1230-1231
###################################################
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, subset = SUBSET, method = "exact.network", mcontrol = mControl(tsmethod = "abs"))


###################################################
### code chunk number 31: intervalCensoring.Rnw:1234-1241
###################################################
SUBSET2<-c(1:12, 47:58)
system.time(
independence_test(Surv(left, right, type = "interval2")~treatment, data = bcos, subset = SUBSET2, ytrafo = wlr_trafo, distribution = exact())
)
system.time(
ictest(Surv(left, right, type = "interval2")~treatment, data = bcos, subset = SUBSET2, method = "exact.network", mcontrol = mControl(tsmethod = "abs"))
)


###################################################
### code chunk number 32: intervalCensoring.Rnw:1253-1258
###################################################
L<-c(2, 5, 1, 1, 9, 8, 10)
R<-c(3, 6, 7, 7, 12, 10, 13)
group<-c(0, 0, 1, 1, 0, 1, 0)
example1<-data.frame(L, R, group)
example1


###################################################
### code chunk number 33: intervalCensoring.Rnw:1278-1279
###################################################
summary(icfit(L, R), digits = 12)


###################################################
### code chunk number 34: intervalCensoring.Rnw:1282-1283
###################################################
print(3/14, digits = 12)


###################################################
### code chunk number 35: intervalCensoring.Rnw:1308-1312
###################################################
score1<-wlr_trafo(Surv(L, R, type = "interval2"))
cm<-chooseMatrix(7, 3)
T<- ( (1-cm) %*% score1 )/4 - ( cm %*% score1 )/3   
cbind(cm, T)[order(T), ][1:9, ]


