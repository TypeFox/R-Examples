### R code from vignette source 'partDSAvig.Rnw'

###################################################
### code chunk number 1: partDSAvig.Rnw:144-147
###################################################
y.out=as.factor(sample(c("a", "b", "c"), 50, TRUE) )
x1=rexp(50)
x2=runif(50)


###################################################
### code chunk number 2: partDSAvig.Rnw:153-155
###################################################
library(partDSA)
#model1<-partDSA(x=data.frame(x1,x2),y=y.out)


###################################################
### code chunk number 3: partDSAvig.Rnw:163-166
###################################################
y.out.test=as.factor(sample(c("a", "b", "c"), 100, TRUE) )
x1.test=rexp(100)
x2.test=runif(100)


###################################################
### code chunk number 4: partDSAvig.Rnw:170-171
###################################################
model2<-partDSA(x=data.frame(x1,x2),y=y.out,x.test=data.frame(x1=x1.test,x2=x2.test),y.test=y.out.test)


###################################################
### code chunk number 5: partDSAvig.Rnw:358-359
###################################################
model4<-partDSA(x=data.frame(x1,x2),y=y.out,control=DSA.control(missing="no",cut.off.growth=2))


###################################################
### code chunk number 6: partDSAvig.Rnw:365-366
###################################################
print(model4)


###################################################
### code chunk number 7: partDSAvig.Rnw:459-462
###################################################
data("GBSG2", package = "TH.data")
mdl1<-partDSA(x=data.frame(GBSG2[,c(1:8)]),y=as.factor(GBSG2$cens),control=DSA.control(cut.off.growth=5,loss.function="gini",minsplit=30,minbuck=10))
print(mdl1)


