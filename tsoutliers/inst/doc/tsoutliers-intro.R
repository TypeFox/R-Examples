### R code from vignette source 'tsoutliers-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tsoutliers-intro.Rnw:57-62
###################################################
library("tsoutliers")
resNile1 <- tso(y = Nile, types = c("AO", "LS", "TC"),
  tsmethod = "stsm", args.tsmodel = list(model = "local-level"))
resNile1$fit$call$xreg <- NULL
resNile1


###################################################
### code chunk number 2: tsoutliers-intro.Rnw:69-73
###################################################
resNile2 <- tso(y = Nile, types = c("AO", "LS", "TC"),
  remove.method = "bottom-up", tsmethod = "auto.arima",
  args.tsmethod = list(allowdrift = FALSE, ic = "bic"))
resNile2


