### R code from vignette source 'exact2x2.Rnw'

###################################################
### code chunk number 1: exact2x2.Rnw:23-24
###################################################
library(exact2x2)


###################################################
### code chunk number 2: exact2x2.Rnw:32-34
###################################################
 x<-matrix(c(6,12,12,5),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))
x


###################################################
### code chunk number 3: exact2x2.Rnw:37-39
###################################################
ft<-fisher.test(x)
ft


###################################################
### code chunk number 4: exact2x2.Rnw:51-53
###################################################
exact2x2(x,tsmethod="minlike")
exact2x2(x,tsmethod="central")


###################################################
### code chunk number 5: exact2x2.Rnw:59-60
###################################################
exact2x2(x,tsmethod="blaker")


