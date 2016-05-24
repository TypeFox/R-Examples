### R code from vignette source 'regsel-heartdisease3.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease3.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: regsel-heartdisease3.Rnw:17-18 (eval = FALSE)
###################################################
## library(lqa)


###################################################
### code chunk number 3: regsel-heartdisease3.Rnw:21-44 (eval = FALSE)
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
### code chunk number 4: regsel-heartdisease3.Rnw:52-62 (eval = FALSE)
###################################################
## lambda2 <- 0.05
## 
## ### COEF BUILD-UPS
## 
## main.text<-"Lasso Improved Correlation Based"
## penalty.family<-licb
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## offset.values = c (NA, lambda2),add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,
## show.standardized=TRUE)


###################################################
### code chunk number 5: regsel-heartdisease3.Rnw:66-73 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.lab,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)
## 


###################################################
### code chunk number 6: regsel-heartdisease3.Rnw:81-86 (eval = FALSE)
###################################################
## main.text<-"Correlation Based"
## penalty.family<-penalreg
## 
## Plot.mat<-plot.lqa (y = y, x = X, family=family, penalty.family=penalty.family,
## add.MLE = FALSE, ret.true=TRUE,really.plot = FALSE,show.standardized=TRUE,gamma=0.01)


###################################################
### code chunk number 7: regsel-heartdisease3.Rnw:90-97 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(Plot.mat$s1,Plot.mat$beta.mat,type="l",ylab=ylab.text,xlab=xlab.text,
## main=main.text,lwd=size.width)
## axis(4, at = Plot.mat$beta.mat[1, ], labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)
## 


###################################################
### code chunk number 8: regsel-heartdisease3.Rnw:102-107 (eval = FALSE)
###################################################
## parcorr=10.5
## max.steps=20
## Path<- ForwardBoost (X.std, y,family = binomial(), penalty = penalreg(parcorr),
## intercept =  TRUE,   nu = 1, monotonic = TRUE,control=lqa.control(max.steps=20,
## conv.stop = FALSE))


###################################################
### code chunk number 9: regsel-heartdisease3.Rnw:110-117 (eval = FALSE)
###################################################
## 
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path$beta.mat[,-1]))/max(rowSums(abs(Path$beta.mat[,-1]))),
## Path$beta.mat[,-1]*sqrt(n),type="l",ylab=ylab.text,xlab=xlab.text,main=
##   "ForwardBoost",lwd=2.0)
## axis(4, at = Path$beta.mat[max.steps, -1]*sqrt(n), labels = colnames(X),
## adj = 0, las = 1,cex.axis=size.right)


