### R code from vignette source 'stsm-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: stsm-intro.Rnw:42-46
###################################################
library("stsm")
m <- stsm.model(model = "BSM", y = log(AirPassengers), transPars = "StructTS")
res <- maxlik.td.optim(m = m, KF.args = list(P0cov = TRUE), method = "L-BFGS-B")
round(get.pars(res$model) * 10000, 3)


