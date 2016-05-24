### R code from vignette source 'strucchange-intro.Rnw'

###################################################
### code chunk number 1: data
###################################################
library("strucchange")
data("USIncExp")
plot(USIncExp, plot.type = "single", col = 1:2, ylab = "billion US$")
legend(1960, max(USIncExp), c("income", "expenditures"),
       lty = c(1,1), col = 1:2, bty = "n")


###################################################
### code chunk number 2: subset
###################################################
library("strucchange")
data("USIncExp")
USIncExp2 <- window(USIncExp, start = c(1985,12))


###################################################
### code chunk number 3: ecm-setup
###################################################
coint.res <- residuals(lm(expenditure ~ income, data = USIncExp2))
coint.res <- lag(ts(coint.res, start = c(1985,12), freq = 12), k = -1)
USIncExp2 <- cbind(USIncExp2, diff(USIncExp2), coint.res)
USIncExp2 <- window(USIncExp2, start = c(1986,1), end = c(2001,2))
colnames(USIncExp2) <- c("income", "expenditure", "diff.income",
                         "diff.expenditure", "coint.res")
ecm.model <- diff.expenditure ~ coint.res + diff.income


###################################################
### code chunk number 4: ts-used
###################################################
plot(USIncExp2[,3:5], main = "")


###################################################
### code chunk number 5: efp
###################################################
ocus <- efp(ecm.model, type="OLS-CUSUM", data=USIncExp2)
me <- efp(ecm.model, type="ME", data=USIncExp2, h=0.2)


###################################################
### code chunk number 6: efp-boundary
###################################################
bound.ocus <- boundary(ocus, alpha=0.05)


###################################################
### code chunk number 7: OLS-CUSUM
###################################################
plot(ocus)


###################################################
### code chunk number 8: efp-boundary2 (eval = FALSE)
###################################################
## plot(ocus, boundary = FALSE)
## lines(bound.ocus, col = 4)
## lines(-bound.ocus, col = 4)


###################################################
### code chunk number 9: ME-null
###################################################
plot(me, functional = NULL)


###################################################
### code chunk number 10: efp-sctest (eval = FALSE)
###################################################
## sctest(ocus)


###################################################
### code chunk number 11: efp-sctest2
###################################################
sctest(ecm.model, type="OLS-CUSUM", data=USIncExp2)


###################################################
### code chunk number 12: Fstats
###################################################
fs <- Fstats(ecm.model, from = c(1990, 1), to = c(1999,6), data = USIncExp2)


###################################################
### code chunk number 13: Fstats-plot
###################################################
plot(fs)


###################################################
### code chunk number 14: pval-plot (eval = FALSE)
###################################################
## plot(fs, pval=TRUE)


###################################################
### code chunk number 15: aveF-plot (eval = FALSE)
###################################################
## plot(fs, aveF=TRUE)


###################################################
### code chunk number 16: Fstats-sctest (eval = FALSE)
###################################################
## sctest(fs, type="expF")


###################################################
### code chunk number 17: Fstats-sctest2
###################################################
sctest(ecm.model, type = "expF", from = 49, to = 162, data = USIncExp2)


###################################################
### code chunk number 18: mefp
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.mefp <- mefp(ecm.model, type = "ME", data = USIncExp3, alpha = 0.05)


###################################################
### code chunk number 19: monitor1
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1990,12))
me.mefp <- monitor(me.mefp)


###################################################
### code chunk number 20: monitor2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)
me.mefp


###################################################
### code chunk number 21: monitor-plot
###################################################
plot(me.mefp)


###################################################
### code chunk number 22: mefp2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.efp <- efp(ecm.model, type = "ME", data = USIncExp3, h = 0.5)
me.mefp <- mefp(me.efp, alpha=0.05)


###################################################
### code chunk number 23: monitor3
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)


###################################################
### code chunk number 24: monitor-plot2
###################################################
plot(me.mefp)


