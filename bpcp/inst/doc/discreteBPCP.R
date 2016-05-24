### R code from vignette source 'discreteBPCP.Rnw'

###################################################
### code chunk number 1: discreteBPCP.Rnw:50-51
###################################################
library(bpcp)


###################################################
### code chunk number 2: discreteBPCP.Rnw:113-114
###################################################
 qbeta(c(.025,.975),1,4)


###################################################
### code chunk number 3: discreteBPCP.Rnw:129-133
###################################################
library(bpcp)
packageVersion("bpcp")
b2<-bpcp(c(3,7,8,14),c(1,1,1,1))
summary(b2)


