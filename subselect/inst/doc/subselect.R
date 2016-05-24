### R code from vignette source 'subselect.Rnw'

###################################################
### code chunk number 1: subselect.Rnw:182-183
###################################################
library(subselect) 


###################################################
### code chunk number 2: subselect.Rnw:348-350
###################################################
rm.coef(mat=var(iris[,-5]),indices=c(3,4))



###################################################
### code chunk number 3: subselect.Rnw:353-355
###################################################
rm.coef(var(iris[,-5]),c(3,4))^2



###################################################
### code chunk number 4: subselect.Rnw:369-370
###################################################
rm.coef(var(iris[,-5]), indices=matrix(nrow=2,ncol=3,byrow=TRUE,c(1,2,3,1,2,4)))


###################################################
### code chunk number 5: subselect.Rnw:384-390
###################################################
subsets  <- array(data=c(3,2,0,0,0,0,1,1,2,2,3,4), dim=c(2,3,2))
colnames(subsets) <- paste("V",1:3,sep="")
rownames(subsets) <- paste("Solution",1:2)
dimnames(subsets)[[3]]<-paste("Size",c(1,3))
subsets
rm.coef(var(iris[,-5]),indices=subsets)


###################################################
### code chunk number 6: subselect.Rnw:478-480
###################################################
gcd.coef(var(iris[,-5]),ind=c(3,4),pcind=c(1,2))



###################################################
### code chunk number 7: subselect.Rnw:485-487
###################################################
gcd.coef(var(iris[,-5]),ind=c(1,2))



###################################################
### code chunk number 8: subselect.Rnw:553-555
###################################################
data(farm)
rv.coef(cor(farm),ind=c(2,37,57,59))


###################################################
### code chunk number 9: subselect.Rnw:563-565
###################################################
rv.coef(cor(farm), indices=matrix(nrow=2,ncol=4,byrow=TRUE,c(2,12,56,59,2,3,11,59)))



###################################################
### code chunk number 10: subselect.Rnw:819-820
###################################################
lmHmat(x=iris[,2:4], y=iris[,1])


###################################################
### code chunk number 11: subselect.Rnw:862-863
###################################################
ldaHmat(x=iris[,1:4], grouping=iris$Species)


###################################################
### code chunk number 12: subselect.Rnw:869-872
###################################################
attach(iris)
ldaHmat(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width)
detach(iris)


###################################################
### code chunk number 13: subselect.Rnw:932-940
###################################################
library(MASS)
data(crabs)
lFL <- log(crabs$FL) ; lRW <- log(crabs$RW); lCL <- log(crabs$CL); lCW <- log(crabs$CW)
C <- matrix(0.,nrow=2,ncol=4)
C[1,3] = C[2,4] = 1.
C
Hmat5 <- glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp*sex,C=C,data=crabs)
Hmat5


###################################################
### code chunk number 14: subselect.Rnw:985-989
###################################################
library(MASS)
data(Cars93)
CarsHmat <- lmHmat(x=Cars93[c(7:8,12:15,17:22,25)],y=Cars93[5])
ccr12.coef(mat=CarsHmat$mat, H=CarsHmat$H, r=CarsHmat$r, indices=c(4,5,10,11))


###################################################
### code chunk number 15: subselect.Rnw:1032-1034
###################################################
irisHmat <- ldaHmat(iris[1:4],iris$Species)
tau2.coef(irisHmat$mat,H=irisHmat$H,r=irisHmat$r,c(1,3))


###################################################
### code chunk number 16: subselect.Rnw:1073-1075
###################################################
irisHmat <- ldaHmat(iris[1:4],iris$Species)
xi2.coef(irisHmat$mat,H=irisHmat$H,r=irisHmat$r,c(1,3))


###################################################
### code chunk number 17: subselect.Rnw:1113-1115
###################################################
irisHmat <- ldaHmat(iris[1:4],iris$Species)
zeta2.coef(irisHmat$mat,H=irisHmat$H,r=irisHmat$r,c(1,3))


###################################################
### code chunk number 18: subselect.Rnw:1183-1188
###################################################
iris2sp <- iris[iris$Species != "setosa",]
modelfit <- glm(Species ~ Sepal.Length + Sepal.Width + Petal.Length +
Petal.Width, data=iris2sp, family=binomial)  
Hmat <- glmHmat(modelfit)
Hmat


###################################################
### code chunk number 19: subselect.Rnw:1231-1239
###################################################
library(MASS) 
lFL <- log(crabs$FL)
lRW <- log(crabs$RW)
lCL <- log(crabs$CL)
lCW <- log(crabs$CW)
logrfit <- glm(sex ~ FL + RW + CL + CW  + lFL + lRW + lCL + lCW,data=crabs,family=binomial)
lHmat <- glmHmat(logrfit) 
wald.coef(lHmat$mat,lHmat$H,indices=c(1,6,7),tolsym=1E-06)  


###################################################
### code chunk number 20: subselect.Rnw:1406-1408
###################################################
data(swiss)
eleaps(cor(swiss),nsol=3, criterion="RM")


###################################################
### code chunk number 21: subselect.Rnw:1426-1429
###################################################
data(swiss)
swiss.gcd <- eleaps(cor(swiss),kmin=2,kmax=3,exclude=6,include=1,nsol=3,criterion="gcd",pcindices=1:3)
swiss.gcd


###################################################
### code chunk number 22: subselect.Rnw:1436-1437
###################################################
rm.coef(mat=cor(swiss), indices=swiss.gcd$subsets)


###################################################
### code chunk number 23: subselect.Rnw:1448-1450
###################################################
irisHmat <- ldaHmat(iris[1:4],iris$Species)
eleaps(irisHmat$mat,kmin=2,kmax=3,H=irisHmat$H,r=irisHmat$r,crit="ccr12")


###################################################
### code chunk number 24: subselect.Rnw:1471-1481
###################################################
library(MASS)
data(Cars93)
Cars93.xgroup <- Cars93[,c(7:8,12:15,17:22,25)]
CarsHmat <- lmHmat(Cars93.xgroup,Cars93[,c(4,6)])
colnames(Cars93[,c(4,6)])
colnames(Cars93.xgroup)
#colnames(CarsHmat$mat)
Cars.eleaps <- eleaps(CarsHmat$mat, kmin=4, kmax=6, H=CarsHmat$H, r=CarsHmat$r, crit="zeta2", tolsym=1e-9)
Cars.eleaps$bestvalues
Cars.eleaps$bestsets


###################################################
### code chunk number 25: subselect.Rnw:1492-1496
###################################################
iris2sp <- iris[iris$Species != "setosa",]
logrfit <- glm(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,iris2sp,family=binomial)
Hmat <- glmHmat(logrfit)
eleaps(Hmat$mat, H=Hmat$H, r=Hmat$r, criterion="Wald", nsol=3)


###################################################
### code chunk number 26: subselect.Rnw:1606-1608
###################################################
data(swiss)
anneal(cor(swiss),kmin=2,kmax=3,nsol=4,niter=10,criterion="RM")


###################################################
### code chunk number 27: subselect.Rnw:1613-1615
###################################################
data(farm)
anneal(cor(farm), kmin=6, nsol=5, criterion="rv")


###################################################
### code chunk number 28: subselect.Rnw:1632-1638
###################################################
library(ISwR)
cystfibrHmat <- lmHmat(pemax ~ age+sex+height+weight+bmp+fev1+rv+frc+tlc, data=cystfibr) 
colnames(cystfibrHmat$mat)
cystfibr.tau2 <- anneal(cystfibrHmat$mat, kmin=4, kmax=6, H=cystfibrHmat$H, r=cystfibrHmat$r, crit="tau2")  
cystfibr.tau2$bestvalues  
cystfibr.tau2$bestsets  


###################################################
### code chunk number 29: subselect.Rnw:1651-1652
###################################################
summary(lm(pemax ~ weight+bmp+fev1+rv, data=cystfibr))$r.squared


###################################################
### code chunk number 30: subselect.Rnw:1661-1662
###################################################
xi2.coef(mat=cystfibrHmat$mat, indices=cystfibr.tau2$bestsets,  H=cystfibrHmat$H, r=cystfibrHmat$r)


###################################################
### code chunk number 31: subselect.Rnw:1674-1681
###################################################
data(iris)
iris2sp <- iris[iris$Species != "setosa",]
logrfit <- glm(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,iris2sp,family=binomial)
Hmat <- glmHmat(logrfit)
iris2p.Wald <- anneal(Hmat$mat,1,3,H=Hmat$H,r=1,nsol=5,criterion="Wald")
iris2p.Wald$bestsets
iris2p.Wald$bestvalues


###################################################
### code chunk number 32: subselect.Rnw:1811-1814
###################################################
farm.gcd <- genetic(cor(farm), kmin=10, crit="gcd")
farm.gcd$bestsets
farm.gcd$bestvalues


###################################################
### code chunk number 33: subselect.Rnw:1822-1823
###################################################
unique(farm.gcd$subsets)


###################################################
### code chunk number 34: subselect.Rnw:1835-1837
###################################################
dim(unique(genetic(cor(farm), kmin=10, crit="gcd")$subsets))
dim(unique(genetic(cor(farm), kmin=10, maxclone=0, crit="gcd")$subsets))


###################################################
### code chunk number 35: subselect.Rnw:1860-1867
###################################################
data(farm) 
farm.xgroup <- farm[,-c(1,2,3,4)]
farmHmat <- lmHmat(farm.xgroup,farm[,1:4])
colnames(farmHmat$mat)
farm.gen <- genetic(farmHmat$mat, kmin=4, kmax=6, H=farmHmat$H, r=farmHmat$r,crit="zeta2", maxclone=0, popsize=150)
farm.gen$bestvalues
farm.gen$bestsets


###################################################
### code chunk number 36: subselect.Rnw:1974-1980
###################################################
swiss.imp1 <- improve(mat=cor(swiss),kmin=2,kmax=3,nsol=4,criterion="GCD")
swiss.imp2 <- improve(cor(swiss),2,3,nsol=4,criterion="GCD",include=c(1),exclude=6)
swiss.imp1$bestvalues
swiss.imp1$bestsets
swiss.imp2$bestvalues
swiss.imp2$bestsets


###################################################
### code chunk number 37: subselect.Rnw:1992-1995
###################################################
data(iris)
irisHmat <- ldaHmat(iris[1:4],iris$Species)
improve(irisHmat$mat,kmin=2,kmax=3,H=irisHmat$H,r=irisHmat$r,crit="ccr12")


###################################################
### code chunk number 38: subselect.Rnw:2010-2018
###################################################
library(MASS)
data(crabs)
lFL <- log(crabs$FL) ; lRW <- log(crabs$RW); lCL <- log(crabs$CL); lCW <- log(crabs$CW)
C <- matrix(0.,nrow=2,ncol=4)
C[1,3] = C[2,4] = 1.
C
Hmat5 <- glhHmat(cbind(FL,RW,CL,CW,lFL,lRW,lCL,lCW) ~ sp*sex,C=C,data=crabs)
improve(mat=Hmat5$mat, kmin=4, nsol=3, H=Hmat5$H, r=Hmat5$r, crit="xi2",tolsym=1e-06)


