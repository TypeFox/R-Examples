### R code from vignette source 'regsel-heartdisease1.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease1.Rnw:13-15 (eval = FALSE)
###################################################
## options(width=80)
## rm(list=ls(all=TRUE))


###################################################
### code chunk number 2: regsel-heartdisease1.Rnw:19-41 (eval = FALSE)
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
### code chunk number 3: regsel-heartdisease1.Rnw:44-46 (eval = FALSE)
###################################################
## library(lqa)
## library(glmnet)


###################################################
### code chunk number 4: regsel-heartdisease1.Rnw:52-57 (eval = FALSE)
###################################################
## main.text<-"Ridge"
## penalty.family<-ridge
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,show.standardized=TRUE,gamma=0.01)


###################################################
### code chunk number 5: regsel-heartdisease1.Rnw:60-66 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width, col=colour)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 6: regsel-heartdisease1.Rnw:72-77 (eval = FALSE)
###################################################
## main.text<-"Lasso"
## penalty.family<-lasso
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## add.MLE = TRUE, ret.true=TRUE,really.plot = FALSE,show.standardized=TRUE,gamma=0.5)


###################################################
### code chunk number 7: regsel-heartdisease1.Rnw:81-87 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width, col=colour)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 8: regsel-heartdisease1.Rnw:92-107 (eval = FALSE)
###################################################
## #Tuning parameter
## g<-1
## ### Weights
## 
## mle.std<-glm(y~X.std,family=binomial(),control=glm.control(maxit=1000))$coef[-1]
## w<-1/(abs(mle.std)^g)
## 
## ### COEF BUILD-UPS
## 
## main.text<-"Adaptive Lasso"
## penalty.family<-adaptive.lasso
## 
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,show.standardized=TRUE,weigths=w)


###################################################
### code chunk number 9: regsel-heartdisease1.Rnw:110-115 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width, col=colour)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 10: regsel-heartdisease1.Rnw:118-123 (eval = FALSE)
###################################################
## opt.min<-Inf
## 
## ##### glmnet for fixed choice
## alpha<-0.5
## Path<-glmnet(X.std, y,family="binomial", alpha = alpha, nlambda = 100)


###################################################
### code chunk number 11: regsel-heartdisease1.Rnw:126-132 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(colSums(abs(Path$beta))/max(colSums(abs(Path$beta))),t(Path$beta)*sqrt(n)
## ,type="l",ylab=ylab.text,xlab=xlab.text,main="Elastic Net (glmnet, a=0.5)",
## lwd=size.width, col=1)
## axis(4, at = Path$beta[,ncol(Path$beta)]*sqrt(n), labels = colnames(X), adj = 0,
##  las = 1,cex.axis=size.right)


###################################################
### code chunk number 12: regsel-heartdisease1.Rnw:135-137 (eval = FALSE)
###################################################
## alpha<-0.2
## Path<-glmnet(X.std, y,family="binomial", alpha = alpha, nlambda = 100)


###################################################
### code chunk number 13: regsel-heartdisease1.Rnw:140-148 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(colSums(abs(Path$beta))/max(colSums(abs(Path$beta))),t(Path$beta)*sqrt(n)
## ,type="l",ylab=ylab.text,xlab=xlab.text,main="Elastic Net (glmnet, a=0.2)",
## lwd=size.width)
## axis(4, at = Path$beta[,ncol(Path$beta)]*sqrt(n), labels = colnames(X), adj = 0,
## las = 1,cex.axis=size.right)
## 


###################################################
### code chunk number 14: regsel-heartdisease1.Rnw:155-164 (eval = FALSE)
###################################################
## lambda2 <- 1
## 
## ### COEF BUILD-UPS
## 
## main.text<-"Elastic Net"
## penalty.family<-enet
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## offset.values = c (NA, lambda2),add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,
## standardize = TRUE,show.standardized=TRUE)


###################################################
### code chunk number 15: regsel-heartdisease1.Rnw:168-173 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


