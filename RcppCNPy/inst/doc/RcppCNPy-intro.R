### R code from vignette source 'RcppCNPy-intro.Rnw'

###################################################
### code chunk number 1: RcppCNPy-intro.Rnw:44-46
###################################################
prettyVersion <- packageDescription("RcppCNPy")[["Version"]]
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 3: RcppCNPy-intro.Rnw:151-157
###################################################
library(methods)
library(RcppCNPy)
M <- matrix(0:11, 3, 4, byrow=TRUE) * 1.1
v <- v <- 0:4 * 1.1
npySave("fmat.npy", M)
npySave("fvec.npy", v)


###################################################
### code chunk number 4: RcppCNPy-intro.Rnw:160-167
###################################################
library(RcppCNPy)

mat <- npyLoad("fmat.npy")
mat

vec <- npyLoad("fvec.npy")
vec


###################################################
### code chunk number 5: RcppCNPy-intro.Rnw:178-179 (eval = FALSE)
###################################################
## mat2 <- npyLoad("fmat.npy.gz")


###################################################
### code chunk number 6: RcppCNPy-intro.Rnw:189-196
###################################################
set.seed(42)
m <- matrix(sort(rnorm(6)), 3, 2)
m
npySave("randmat.npy", m)
v <- seq(10, 12)
v
npySave("simplevec.npy", v)


###################################################
### code chunk number 8: RcppCNPy-intro.Rnw:296-300
###################################################
unlink("fmat.npy")
unlink("fvec.npy")
unlink("randmat.npy")
unlink("simplevec.npy")


