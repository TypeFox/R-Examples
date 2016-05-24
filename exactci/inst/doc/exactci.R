### R code from vignette source 'exactci.Rnw'

###################################################
### code chunk number 1: exactci.Rnw:23-24
###################################################
library(exactci)


###################################################
### code chunk number 2: exactci.Rnw:36-37
###################################################
poisson.test(c(2,10),c(17877,20000))


###################################################
### code chunk number 3: exactci.Rnw:61-65
###################################################
library(exactci)
poisson.exact(c(2,10),c(17877,20000),tsmethod="central")
poisson.exact(c(2,10),c(17877,20000),tsmethod="minlike")
poisson.exact(c(2,10),c(17877,20000),tsmethod="blaker")


