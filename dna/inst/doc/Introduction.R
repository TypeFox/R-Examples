### R code from vignette source 'Introduction.Rnw'

###################################################
### code chunk number 1: Introduction.Rnw:165-168
###################################################
library("dna")
data("HeavyMice")
data("LeanMice")


###################################################
### code chunk number 2: Introduction.Rnw:419-422
###################################################
set.seed(26)
s=matrix(runif(100,-1,1),10,10);diag(s)=1;s=round((s+t(s))/2,1)
s


###################################################
### code chunk number 3: Introduction.Rnw:460-461
###################################################
network.modules(s,m=3,epsilon=.7,plot=TRUE,interactive=FALSE)


