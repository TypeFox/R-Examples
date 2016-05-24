### R code from vignette source 'CNY.Rnw'

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
### code chunk number 3: cny-data
###################################################
cny <- fxreturns("CNY", frequency = "daily",
  start = as.Date("2005-07-25"), end = as.Date("2009-07-31"),
  other = c("USD", "JPY", "EUR", "GBP"), data = FXRatesCHF)


###################################################
### code chunk number 4: cny-fitting
###################################################
cny_lm <- fxlm(CNY ~ USD + JPY + EUR + GBP,
  data = window(cny, end = as.Date("2005-10-31")))
summary(cny_lm)


###################################################
### code chunk number 5: cny-testing
###################################################
cny_efp <- gefp(cny_lm, fit = NULL)


###################################################
### code chunk number 6: cny-hprocess (eval = FALSE)
###################################################
## plot(cny_efp, aggregate = FALSE, ylim = c(-1.85, 1.85))


###################################################
### code chunk number 7: cny-hprocess1
###################################################
plot(cny_efp, aggregate = FALSE, ylim = c(-1.85, 1.85))


###################################################
### code chunk number 8: cny-sctest
###################################################
sctest(cny_efp)


###################################################
### code chunk number 9: cny-monitoring (eval = FALSE)
###################################################
## cny_mon <- fxmonitor(CNY ~ USD + JPY + EUR + GBP,
##   data = window(cny, end = as.Date("2006-05-31")),
##   start = as.Date("2005-11-01"), end = 4)
## plot(cny_mon, aggregate = FALSE)


###################################################
### code chunk number 10: cny-mprocess
###################################################
cny_mon <- fxmonitor(CNY ~ USD + JPY + EUR + GBP,
  data = window(cny, end = as.Date("2006-05-31")),
  start = as.Date("2005-11-01"), end = 4)
plot(cny_mon, aggregate = FALSE)


###################################################
### code chunk number 11: cny-monitorbreak
###################################################
cny_mon


###################################################
### code chunk number 12: cny-dating (eval = FALSE)
###################################################
## cny_reg <- fxregimes(CNY ~ USD + JPY + EUR + GBP,
##   data = cny, h = 20, breaks = 10)


###################################################
### code chunk number 13: cny-dating1
###################################################
if(file.exists("cny_reg.rda")) load("cny_reg.rda") else {
cny_reg <- fxregimes(CNY ~ USD + JPY + EUR + GBP,
  data = cny, h = 20, breaks = 10)
save(cny_reg, file = "cny_reg.rda")
}
if(cleanup) file.remove("cny_reg.rda")


###################################################
### code chunk number 14: cny-breaks (eval = FALSE)
###################################################
## plot(cny_reg)


###################################################
### code chunk number 15: cny-breaks1
###################################################
plot(cny_reg)


###################################################
### code chunk number 16: cny-confint
###################################################
confint(cny_reg, level = 0.9)


###################################################
### code chunk number 17: cny-coef
###################################################
coef(cny_reg)


###################################################
### code chunk number 18: cny-refit
###################################################
cny_rf <- refit(cny_reg)
lapply(cny_rf, summary)


