### R code from vignette source 'rtkpp-introduction.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(rtkpp)
rtkpp.version <- packageDescription("rtkpp")$Version
rtkpp.date <- packageDescription("rtkpp")$Date


###################################################
### code chunk number 2: rtkpp-introduction.Rnw:82-83
###################################################
.Call("stk_version", FALSE, PACKAGE="rtkpp")


