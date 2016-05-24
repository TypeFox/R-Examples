### R code from vignette source 'overview.Rnw'

###################################################
### code chunk number 1: a1
###################################################
library(dr)
opt <- options(width=66)


###################################################
### code chunk number 2: a11
###################################################
summary(s0 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),data=ais,slice.function=dr.slices.arc,nslices=8,
  chi2approx="wood",numdir=4,method="sir"))


###################################################
### code chunk number 3: a2
###################################################
dr.coordinate.test(s0,hypothesis=~.-log(RCC))


###################################################
### code chunk number 4: a2
###################################################
dr.coordinate.test(s0,hypothesis=~.-log(RCC),d=2)


###################################################
### code chunk number 5: a2
###################################################
m0 <- drop1(s0,update=TRUE)


###################################################
### code chunk number 6: overview.Rnw:563-564
###################################################
s1a <- dr.step(s0,scope=~log(Wt),stop=0.20)


###################################################
### code chunk number 7: overview.Rnw:599-600
###################################################
summary(s1 <- update(s0, group=~Sex))


###################################################
### code chunk number 8: a2
###################################################
s2 <- update(s0,method="save")
summary(s2)


###################################################
### code chunk number 9: save
###################################################
drop1(s1,update=FALSE)


###################################################
### code chunk number 10: overview.Rnw:706-707
###################################################
summary(s3 <- update(s2,group=~Sex))


###################################################
### code chunk number 11: overview.Rnw:745-746
###################################################
summary(s2 <- update(s0,method="phdres"))


###################################################
### code chunk number 12: one
###################################################
(m1 <- dr(LBM~log(Ht)+log(Wt)+log(SSF)+log(RCC)+log(WCC)+log(Ferr)+
           log(Hc)+log(Hg),data=ais,method="ire",nslices=8,numdir=4,
           slice.function=dr.slices.arc,itmax=200,steps=1,eps=1.e-6))


###################################################
### code chunk number 13: two
###################################################
dr.basis(m1,numdir=2)


###################################################
### code chunk number 14: three
###################################################
dr.basis(m1,3)


###################################################
### code chunk number 15: mct
###################################################
dr.coordinate.test(m1,~.-log(Hg))
dr.coordinate.test(m1,~.-log(Hg),d=2)


###################################################
### code chunk number 16: drop1
###################################################
drop1(m1,update=FALSE)


###################################################
### code chunk number 17: pire
###################################################
m2 <- dr(LBM~log(Ht)+log(Wt)+log(SSF)+log(RCC)+log(WCC)+log(Ferr)+
           log(Hc)+log(Hg),group=~Sex,data=ais,method="ire",nslices=8,
           numdir=4,slice.function=dr.slices.arc,itmax=200,steps=1,
           eps=1.e-6)


###################################################
### code chunk number 18: pire1
###################################################
m2


###################################################
### code chunk number 19: overview.Rnw:1102-1104
###################################################
wts <- dr.weights(LBM~Ht+Wt+RCC+WCC,data=ais)
i1 <- dr(LBM~Ht+Wt+RCC+WCC,weights=wts,method="phdres",data=ais)


###################################################
### code chunk number 20: overview.Rnw:1113-1115
###################################################
y1 <- c(1,1,1,2,3,4,5,6,7,8,8,8)
dr.slices(y1,3)


###################################################
### code chunk number 21: overview.Rnw:1128-1130
###################################################
y2 <- c(1,2,3,4,1,2,3,4,1,2,3,4)
dr.slices(cbind(y1,y2),5)


###################################################
### code chunk number 22: overview.Rnw:1145-1146
###################################################
dr.permutation.test(s0,npermute=99,numdir=4)


