#15-1-2014 MRC-Epid JHZ

library(gap)
library(kinship2)

ped <- with(ped51,pedigree(id,fid,mid,sex,aff=bt))
plot(ped)
kmat <- with(ped51,kinship(id, fid, mid))
k2 <- as.matrix(2*kmat)
N <- dim(l51)[1]
I <- diag(N)
dimnames(I) <- dimnames(k2)
k3 <- k2+I
detach(package:kinship2)

library(regress)

regress(qt ~ 1, ~I, data=ped51)
regress(qt ~ 1, ~k2, data=ped51)
regress(qt ~ 1, ~k2+I, data=ped51)
regress(qt ~ 1, ~k3, data=ped51)

library(coxme)

lmekin(qt ~ 1 + (1|id), data=ped51)
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k2))
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k2,I))
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k3))

coxme(Surv(qt,bt)~1+(1|id),data=ped51)
coxme(Surv(qt,bt)~1+(1|id),data=ped51,varlist=list(k2))
coxme(Surv(qt,bt)~1+(1|id),data=ped51,varlist=list(k3))

library(pedigreemm)

ped <- with(ped51,pedigree(fid,mid,as.character(id)))
pm0 <- pedigreemm(bt~(1|id), ped51, family="binomial")
pm1 <- pedigreemm(bt~(1|id), ped51, family="binomial", pedigree=list(id=ped))
r1 <- resid(pm1)
l1 <- lm(r1~qt,data=ped51)
summary(l1)
detach(package:pedigreemm)

library(pls)
library(POET)

with(ped51,summary(qt))
m0 <- lm(qt ~ 1, data=ped51)
summary(m0)
m1 <- lm(qt ~ r, data=ped51)
summary(m1)
pca <- prcomp(k2)
screeplot(pca,25)
biplot(pca)
pred <- predict(pca)
l1 <- lm(qt~1+pred[,1:6],data=ped51)
summary(l1)
l2 <- lm(qt~1+pred[,1:35],data=ped51)
summary(l2)
  
p0 <- pcr(qt ~ pred, 25, data=ped51, validation="CV")
summary(p0)
coef(p0)
coefplot(p0)
p1 <- pcr(qt ~ k2, 35, data=ped51, validation="CV")
summary(p1)
coef(p0)
coefplot(p1)
scores(p1)
poet <- POET(k2)
pred <- k2%*%poet$loadings
t0 <- lm(qt~pred,data=ped51)
summary(t0)
poet <- POET(k2,4)
pred <- k2%*%poet$loadings
t1 <- lm(qt~pred,data=ped51)
summary(t1)
