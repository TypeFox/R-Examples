### R code from vignette source 'ZIM.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70)


###################################################
### code chunk number 2: zim
###################################################
library(ZIM)

data(syph)
count <- syph$a33
ar1 <- bshift(count > 0)
trend <- 1:length(count) / 1000


###################################################
### code chunk number 3: zip
###################################################
zim(count ~ ar1 + trend | trend)


###################################################
### code chunk number 4: zinb
###################################################
zim(count ~ ar1 + trend | trend, dist = "zinb")


