### R code from vignette source 'Modalclust.Rnw'

###################################################
### code chunk number 1: Modalclust.Rnw:72-75
###################################################
library(Modalclust)
data(disc2d)
plot(disc2d)


###################################################
### code chunk number 2: Modalclust.Rnw:93-94
###################################################
disc2d.hmac=phmac(disc2d,npart=1,parallel=F)


###################################################
### code chunk number 3: Modalclust.Rnw:112-114
###################################################
disc2d.phmac=phmac(disc2d,npart=2,parallel=T) 
          #npart can be made equal to the number of processors the computer has.


###################################################
### code chunk number 4: summary
###################################################
summary(disc2d.phmac)


###################################################
### code chunk number 5: Modalclust.Rnw:146-147
###################################################
plot.hmac(disc2d.hmac)


###################################################
### code chunk number 6: Modalclust.Rnw:155-160
###################################################
set.seed(20)
mix4=data.frame(rbind(rmvnorm(20,rep(0,4)), rmvnorm(20,rep(2,4)),
                      rmvnorm(20,rep(10,4)),rmvnorm(20,rep(13,4))))
mix4.hmac=phmac(mix4,npart=1)
plot(mix4.hmac)


###################################################
### code chunk number 7: Modalclust.Rnw:166-167
###################################################
plot(mix4.hmac,userclus=rep(c(1,2,3,4),each=20))


###################################################
### code chunk number 8: Modalclust.Rnw:173-174
###################################################
plot(mix4.hmac,n.cluster=2)


###################################################
### code chunk number 9: Modalclust.Rnw:183-187
###################################################
par(mfrow=c(1,2))
hard.hmac(disc2d.hmac,n.cluster=2)
soft.hmac(disc2d.hmac,level=3)
par(mfrow=c(1,1))


###################################################
### code chunk number 10: Modalclust.Rnw:191-195
###################################################
member=hard.hmac(disc2d.hmac,n.cluster=2,plot=F)
member
member.soft=soft.hmac(disc2d.hmac,level=3,plot=F)
head(member.soft$post.prob)


###################################################
### code chunk number 11: Modalclust.Rnw:208-209
###################################################
choose.cluster(disc2d.hmac,n.cluster=2,x=c(0,0))


###################################################
### code chunk number 12: Modalclust.Rnw:217-218
###################################################
contour.hmac(disc2d.hmac,n.cluster=2,col=gray(0.7))


###################################################
### code chunk number 13: Modalclust.Rnw:243-246
###################################################
data(cta20)
data(cta20.hmac)
plot(cta20,pch=16)


###################################################
### code chunk number 14: Modalclust.Rnw:266-271
###################################################
data(logcta20.hmac)
par(mfrow=c(1,2))
hard.hmac(cta20.hmac,n.cluster=3)
hard.hmac(logcta20.hmac,n.cluster=3)
par(mfrow=c(1,1))


###################################################
### code chunk number 15: Modalclust.Rnw:277-283
###################################################
data(oned.hmac)
data(oned)
par(mfrow=c(1,2))
hist(oned,n=20,col="lavender")
hard.hmac(oned.hmac,level=3)
par(mfrow=c(1,1))


###################################################
### code chunk number 16: Modalclust.Rnw:295-297
###################################################
iris.hmac=phmac(iris[,-5],npart=1,parallel=F)
plot(iris.hmac,sep=.02)


###################################################
### code chunk number 17: Modalclust.Rnw:301-302
###################################################
hard.hmac(iris.hmac,n.cluster=2)


