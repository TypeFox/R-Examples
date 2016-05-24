### R code from vignette source 'regsel-heartdisease4.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease4.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: regsel-heartdisease4.Rnw:17-20 (eval = FALSE)
###################################################
## library(lpSolve)
## library(lqa)
## source("gdscode.txt")


###################################################
### code chunk number 3: regsel-heartdisease4.Rnw:23-45 (eval = FALSE)
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
### code chunk number 4: regsel-heartdisease4.Rnw:55-65 (eval = FALSE)
###################################################
## s <- 3
## 
## ### COEF BUILD-UPS
## 
## main.text<-"SCAD"
## penalty.family<-scad
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## offset.values = c (NA, s),add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,
## show.standardized=TRUE,gamma=0.01)


###################################################
### code chunk number 5: regsel-heartdisease4.Rnw:68-75 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width, col=colour)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)
## 


###################################################
### code chunk number 6: regsel-heartdisease4.Rnw:80-86 (eval = FALSE)
###################################################
## Path<-matrix(0,60,p)
## lambda1<-exp(seq(-6, 1, length = 60))
## for(i in 60:1)
## {
## Path[i,]<-dd(y,X.std,lambda=lambda1[i],family=binomial)$gds.beta[-1]*sqrt(n)
## }


###################################################
### code chunk number 7: regsel-heartdisease4.Rnw:90-95 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path,type="l",ylab=ylab.text,
## xlab=xlab.text,main="Dantzig Selector",lwd=size.width, col=colour)
## axis(4, at = Path[1, ], labels = colnames(X), adj = 0, las = 1,cex.axis=size.right)


