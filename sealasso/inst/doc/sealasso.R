### R code from vignette source 'sealasso.Rnw'

###################################################
### code chunk number 1: sealasso.Rnw:42-47
###################################################
library(sealasso)
data(diabetes)      # use the diabetes dataset from "lars" package
x <- diabetes$x
y <- diabetes$y
sealasso(x, y, method = "sealasso")


###################################################
### code chunk number 2: interaction_effect
###################################################
# with quadratic terms
x2 <- cbind(diabetes$x1,diabetes$x2)
object <- sealasso(x2, y)
summary(object)


