### R code from vignette source 'mnps.rnw'

###################################################
### code chunk number 1: mnps.rnw:74-75
###################################################
options(width=60)


###################################################
### code chunk number 2: mnps.rnw:96-99
###################################################
library(twang)
data(AOD)
set.seed(1)


###################################################
### code chunk number 3: mnps.rnw:145-151
###################################################
mnps.AOD <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                 data = AOD, 
                 estimand = "ATE", 
                 verbose = FALSE, 
                 stop.method = c("es.mean", "ks.mean"), 
                 n.trees = 3000)


###################################################
### code chunk number 4: mnps.rnw:202-203
###################################################
    plot(mnps.AOD, plots = 1)


###################################################
### code chunk number 5: mnps.rnw:239-240
###################################################
    plot(mnps.AOD, plots = 2, subset = "es.mean")


###################################################
### code chunk number 6: mnps.rnw:272-273
###################################################
options(width=120)


###################################################
### code chunk number 7: mnps.rnw:277-278
###################################################
    plot(mnps.AOD, plots = 3)


###################################################
### code chunk number 8: mnps.rnw:290-291
###################################################
options(width=120)


###################################################
### code chunk number 9: mnps.rnw:298-299
###################################################
plot(mnps.AOD, plots = 3, pairwiseMax = FALSE, figureRows = 3)


###################################################
### code chunk number 10: mnps.rnw:325-326
###################################################
    plot(mnps.AOD, plots = 4)


###################################################
### code chunk number 11: mnps.rnw:331-332
###################################################
options(width=60)


###################################################
### code chunk number 12: mnps.rnw:352-353
###################################################
    plot(mnps.AOD, plots = 2, subset = "es.mean", singlePlot = 2)


###################################################
### code chunk number 13: mnps.rnw:373-374
###################################################
bal.table(mnps.AOD, digits = 2)


###################################################
### code chunk number 14: mnps.rnw:394-395
###################################################
bal.table(mnps.AOD, collapse.to = 'covariate', digits = 4)


###################################################
### code chunk number 15: mnps.rnw:407-408
###################################################
bal.table(mnps.AOD, collapse.to = 'stop.method', digits = 4)


###################################################
### code chunk number 16: mnps.rnw:436-438
###################################################
bal.table(mnps.AOD, subset.treat = c('community', 'metcbt5'), 
          subset.var = c('white', 'illact', 'crimjust'))


###################################################
### code chunk number 17: mnps.rnw:441-442
###################################################
bal.table(mnps.AOD, subset.stop.method = 'es.mean', collapse.to = 'covariate')


###################################################
### code chunk number 18: mnps.rnw:445-446
###################################################
bal.table(mnps.AOD, es.cutoff = 0.1)


###################################################
### code chunk number 19: mnps.rnw:456-457
###################################################
summary(mnps.AOD)


###################################################
### code chunk number 20: mnps.rnw:461-462
###################################################
options(width=60)


###################################################
### code chunk number 21: mnps.rnw:485-488
###################################################
library(survey)
AOD$w <- get.weights(mnps.AOD, stop.method = "es.mean")
design.mnps <- svydesign(ids=~1, weights=~w, data=AOD)


###################################################
### code chunk number 22: mnps.rnw:493-495
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps)
summary(glm1)


###################################################
### code chunk number 23: mnps.rnw:520-522
###################################################
glm2 <- svyglm(suf12 ~ treat, design = design.mnps, contrast=list(treat=contr.sum))
summary(glm2)


###################################################
### code chunk number 24: mnps.rnw:533-534
###################################################
-sum(coef(glm2)[-1])


###################################################
### code chunk number 25: mnps.rnw:538-539
###################################################
sqrt(c(-1,-1) %*% summary(glm2)$cov.scaled[-1,-1] %*% c(-1,-1))


###################################################
### code chunk number 26: mnps.rnw:561-568
###################################################
mnps.AOD.ATT <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                     data = AOD, 
                     estimand = "ATT", 
                     treatATT = "community", 
                     verbose = FALSE, 
                     n.trees = 3000, 
                     stop.method = c("es.mean", "ks.mean"))


###################################################
### code chunk number 27: mnps.rnw:580-581
###################################################
    plot(mnps.AOD.ATT, plots = 1)


###################################################
### code chunk number 28: mnps.rnw:588-589
###################################################
    plot(mnps.AOD.ATT, plots = 3)


###################################################
### code chunk number 29: mnps.rnw:594-595
###################################################
    plot(mnps.AOD.ATT, plots = 3, pairwiseMax = FALSE)


###################################################
### code chunk number 30: mnps.rnw:603-604
###################################################
    plot(mnps.AOD.ATT, plots = 4)


###################################################
### code chunk number 31: mnps.rnw:608-609
###################################################
options(width=85)


###################################################
### code chunk number 32: mnps.rnw:619-620
###################################################
bal.table(mnps.AOD.ATT, digits = 2)


###################################################
### code chunk number 33: mnps.rnw:623-624
###################################################
bal.table(mnps.AOD.ATT, digits = 2, collapse.to  = "covariate")


###################################################
### code chunk number 34: mnps.rnw:627-628
###################################################
bal.table(mnps.AOD.ATT, digits = 3, collapse.to = "stop.method")


###################################################
### code chunk number 35: mnps.rnw:632-633
###################################################
options(width=60)


###################################################
### code chunk number 36: mnps.rnw:640-643
###################################################
require(survey)
AOD$w.ATT <- get.weights(mnps.AOD.ATT, stop.method = "es.mean")
design.mnps.ATT <- svydesign(ids=~1, weights=~w.ATT, data=AOD)


###################################################
### code chunk number 37: mnps.rnw:647-649
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps.ATT)
summary(glm1)


###################################################
### code chunk number 38: mnps.rnw:743-744
###################################################
options(width=85)


###################################################
### code chunk number 39: mnps.rnw:747-748
###################################################
means.table(mnps.AOD, stop.method = "es.mean", digits = 3)


###################################################
### code chunk number 40: mnps.rnw:759-760
###################################################
means.table(mnps.AOD.ATT, digits = 3)


