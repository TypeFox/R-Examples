### R code from vignette source 'regsel-heartdisease5.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease5.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: regsel-heartdisease5.Rnw:17-19 (eval = FALSE)
###################################################
## library(mboost)
## library(GAMBoost)


###################################################
### code chunk number 3: regsel-heartdisease5.Rnw:22-45 (eval = FALSE)
###################################################
## library(catdata)
## data(heart, package="catdata")
## 
## X<-heart[,-1]
## y<-heart[,1]
## X.std<-scale(X)
## p<-ncol(X)
## n<-length(y)
## family <- binomial()
## n.fold<-10
## 
## 
## ylab.text<-""
## xlab.text<-""
## Width = 6 
## Height = 6 
## oma.vec<-c(1,1,1,3) 
## size.axis=1.4 
## size.lab=1.4 
## size.main=1.4 
## size.right=1.2 
## size.width=2.0
## colour=1


###################################################
### code chunk number 4: regsel-heartdisease5.Rnw:51-54 (eval = FALSE)
###################################################
## Path<-GLMBoost(x=X.std,y,penalty=length(y),standardize=FALSE,family=binomial(),
## stepno=500)
## Path<-Path$beta.linear*sqrt(n)


###################################################
### code chunk number 5: regsel-heartdisease5.Rnw:58-64 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path,type="l",ylab=ylab.text,
## xlab=xlab.text,main="GLMBoost",lwd=size.width, col=colour)
## axis(4,at = Path[500,], labels = colnames(X), adj = 0,las =1,cex.axis=size.right)
## 


###################################################
### code chunk number 6: regsel-heartdisease5.Rnw:68-78 (eval = FALSE)
###################################################
## y.boost<-as.factor(y*2-1)
## X.boost<-X.std/2
## mstop=500
## aux<-glmboost(y.boost~X.boost,family=Binomial(),control=boost_control(mstop=mstop,
## nu=0.1))
## Path<-matrix(0,mstop,p)
## for(i in 1:500)
## {
## Path[i,]<-coef(aux[i],2:10)
## }


###################################################
### code chunk number 7: regsel-heartdisease5.Rnw:82-89 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(c(0,rowSums(abs(Path))/max(rowSums(abs(Path)))),rbind(0,Path*sqrt(n)),
## type="l",ylab=ylab.text,xlab=xlab.text,main="glmboost (mboost)",lwd=size.width,
## col=colour)
## axis(4, at = Path[500,]*sqrt(n), labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


