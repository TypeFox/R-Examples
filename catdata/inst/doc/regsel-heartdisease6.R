### R code from vignette source 'regsel-heartdisease6.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease6.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: regsel-heartdisease6.Rnw:17-18 (eval = FALSE)
###################################################
## library(lqa)


###################################################
### code chunk number 3: regsel-heartdisease6.Rnw:21-46 (eval = FALSE)
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
## 
## library(penalized)


###################################################
### code chunk number 4: regsel-heartdisease6.Rnw:52-63 (eval = FALSE)
###################################################
## alpha <- 0.7
## g<-5
## 
## ### COEF BUILD-UPS
## 
## main.text<-"Generalized Elastic Net"
## penalty.family<-genet
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## offset.values = c (NA, alpha,g),add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,
## show.standardized=TRUE,gamma=0.01)


###################################################
### code chunk number 5: regsel-heartdisease6.Rnw:65-71 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 6: regsel-heartdisease6.Rnw:76-85 (eval = FALSE)
###################################################
## 
## lambda2<-0.1
## 
## main.text<-"Improved Correlation Based"
## penalty.family<-icb
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## offset.values = c (NA, lambda2),add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,
## show.standardized=TRUE,gamma=0.01)


###################################################
### code chunk number 7: regsel-heartdisease6.Rnw:89-95 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 8: regsel-heartdisease6.Rnw:99-108 (eval = FALSE)
###################################################
## Path<-matrix(0,60,p)
## lambda1<-exp(seq(-10, 6, length = 60))
## lambda2 = 1
## 
## for(i in 60:1)
## {
## Path[i,]<-coef(penalized (y,penalized=X.std/sqrt(n), lambda1 = lambda1[i],
## lambda2 = lambda2,model="logistic",standardize=FALSE),"penalized")
## }


###################################################
### code chunk number 9: regsel-heartdisease6.Rnw:111-116 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path,type="l",ylab=ylab.text,
## xlab=xlab.text,main="Enet with penalized",lwd=size.width)
## axis(4, at = Path[1, ], labels = colnames(X), adj = 0, las = 1,cex.axis=size.right)


###################################################
### code chunk number 10: regsel-heartdisease6.Rnw:121-130 (eval = FALSE)
###################################################
## Path<-matrix(0,60,p)
## lambda1<-exp(seq(-10, 6, length = 60))
## lambda2 = 0
## 
## for(i in 60:1)
## {
## Path[i,]<-coef(penalized (y,penalized=X.std/sqrt(n), lambda1 = lambda1[i],
## lambda2 = lambda2,model="logistic",standardize=FALSE),"penalized")
## }


###################################################
### code chunk number 11: regsel-heartdisease6.Rnw:132-136 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path,type="l",ylab=ylab.text,
## xlab=xlab.text,main="Lasso with penalized",lwd=size.width)
## axis(4, at = Path[1, ], labels = colnames(X), adj = 0, las = 1,cex.axis=size.right)


###################################################
### code chunk number 12: regsel-heartdisease6.Rnw:139-140 (eval = FALSE)
###################################################
## detach(package:penalized)


