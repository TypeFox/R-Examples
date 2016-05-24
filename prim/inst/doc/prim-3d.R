### R code from vignette source 'prim-3d.Rnw'

###################################################
### code chunk number 1: prim-3d.Rnw:69-75
###################################################
library(prim)
data(quasiflow)
yflow <- quasiflow[,4]
xflow <- quasiflow[,1:3]
xflowp <- quasiflow[yflow==1,1:3]
xflown <- quasiflow[yflow==-1,1:3]


###################################################
### code chunk number 2: prim-3d.Rnw:79-81
###################################################
pairs(xflowp[1:300,])
pairs(xflown[1:300,])


###################################################
### code chunk number 3: prim-3d.Rnw:88-89
###################################################
pairs(xflowp[1:300,])


###################################################
### code chunk number 4: prim-3d.Rnw:92-93
###################################################
pairs(xflown[1:300,])


###################################################
### code chunk number 5: prim-3d.Rnw:103-105 (eval = FALSE)
###################################################
## qflow.thr <- c(0.38, -0.23)
## qflow.prim <- prim.box(x=xflow, y=yflow, threshold=qflow.thr, threshold.type=0)


###################################################
### code chunk number 6: prim-3d.Rnw:114-115
###################################################
qflow.hdr.pos <- prim.box(x=xflow, y=yflow, threshold=0.38, threshold.type=1)


###################################################
### code chunk number 7: prim-3d.Rnw:118-122
###################################################
qflow.neg <- prim.box(x=xflow, y=yflow, threshold.type=-1)
qflow.hdr.neg1 <- prim.hdr(qflow.neg, threshold=-0.23, threshold.type=-1)
qflow.hdr.neg2 <- prim.hdr(qflow.neg, threshold=-0.43, threshold.type=-1)
qflow.hdr.neg3 <- prim.hdr(qflow.neg, threshold=-0.63, threshold.type=-1)


###################################################
### code chunk number 8: prim-3d.Rnw:126-128
###################################################
qflow.prim2 <- prim.combine(qflow.hdr.pos, qflow.hdr.neg1)
summary(qflow.prim2)


###################################################
### code chunk number 9: prim-3d.Rnw:134-135
###################################################
plot(qflow.prim2)


