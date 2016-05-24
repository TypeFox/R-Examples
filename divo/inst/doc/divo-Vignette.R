### R code from vignette source 'divo-Vignette.Rnw'

###################################################
### code chunk number 1: divo-Vignette.Rnw:26-29
###################################################
library(divo)
data(TCR.Data)
result <- i.inp(x, resample = 25)


###################################################
### code chunk number 2: divo-Vignette.Rnw:33-36
###################################################
data(TCR.Data)
DP <- dp(x, resample = 25)
ENS <- ens(x, resample = 25)


