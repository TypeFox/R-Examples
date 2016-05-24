### R code from vignette source 'expectreg.Rnw'

###################################################
### code chunk number 1: expectreg.Rnw:30-31
###################################################
library(expectreg)


###################################################
### code chunk number 2: expectreg.Rnw:33-35 (eval = FALSE)
###################################################
## help(package = "expectreg")
## data(package = "expectreg")


###################################################
### code chunk number 3: expectreg.Rnw:141-143 (eval = FALSE)
###################################################
## data(india)
## data(dutchboys)


###################################################
### code chunk number 4: expectreg.Rnw:159-160 (eval = FALSE)
###################################################
## data(dutchboys)


###################################################
### code chunk number 5: expectreg.Rnw:162-163 (eval = FALSE)
###################################################
## exp.l <- expectreg.ls(dutchboys[,3] ~ rb(dutchboys[,2],"pspline"),smooth="acv")


###################################################
### code chunk number 6: expectreg.Rnw:177-178 (eval = FALSE)
###################################################
## exp.b <- expectreg.ls(dutchboys[,3] ~ rb(dutchboys[,2],"pspline"),smooth="none",estimate="bundle")


###################################################
### code chunk number 7: expectreg.Rnw:191-192 (eval = FALSE)
###################################################
## exp.r <- expectreg.ls(dutchboys[,3] ~ rb(dutchboys[,2],"pspline"),smooth="schall",estimate="restricted")


###################################################
### code chunk number 8: expectreg.Rnw:205-206 (eval = FALSE)
###################################################
## exp.boost <- expectreg.boost(hgt ~ bbs(age,df=5,degree=2),dutchboys,mstop=rep(500,11)) 


