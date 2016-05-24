### R code from vignette source 'INR.Rnw'

###################################################
### code chunk number 1: cleanup
###################################################
cleanup <- TRUE


###################################################
### code chunk number 2: preliminaries
###################################################
library("fxregime")
data("FXRatesCHF", package = "fxregime")


###################################################
### code chunk number 3: inr-data
###################################################
inr <- fxreturns("INR", frequency = "weekly",
  start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
  other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF)


###################################################
### code chunk number 4: inr-fitting
###################################################
inr_lm <- fxlm(INR ~ USD + JPY + DUR + GBP, data = inr)


###################################################
### code chunk number 5: inr-testing (eval = FALSE)
###################################################
## inr_efp <- gefp(inr_lm, fit = NULL)
## plot(inr_efp, aggregate = FALSE, ylim = c(-1.85, 1.85))


###################################################
### code chunk number 6: inr-hprocess
###################################################
inr_efp <- gefp(inr_lm, fit = NULL)
plot(inr_efp, aggregate = FALSE, ylim = c(-1.85, 1.85))


###################################################
### code chunk number 7: inr-sctest
###################################################
sctest(inr_efp)


###################################################
### code chunk number 8: inr-sctest2
###################################################
sctest(inr_efp, functional = meanL2BB)


###################################################
### code chunk number 9: inr-dating (eval = FALSE)
###################################################
## inr_reg <- fxregimes(INR ~ USD + JPY + DUR + GBP,
##   data = inr, h = 20, breaks = 10)


###################################################
### code chunk number 10: inr-dating1
###################################################
if(file.exists("inr_reg.rda")) load("inr_reg.rda") else {
inr_reg <- fxregimes(INR ~ USD + JPY + DUR + GBP,
  data = inr, h = 20, breaks = 10)
save(inr_reg, file = "inr_reg.rda")
}
if(cleanup) file.remove("inr_reg.rda")


###################################################
### code chunk number 11: inr-breaks (eval = FALSE)
###################################################
## plot(inr_reg)


###################################################
### code chunk number 12: inr-breaks1
###################################################
plot(inr_reg)


###################################################
### code chunk number 13: inr-confint
###################################################
confint(inr_reg, level = 0.9)


###################################################
### code chunk number 14: inr-coef
###################################################
coef(inr_reg)


###################################################
### code chunk number 15: inr-coef
###################################################
inr_rf <- refit(inr_reg)
lapply(inr_rf, summary)


