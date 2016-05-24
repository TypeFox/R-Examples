### R code from vignette source 'regsel-heartdisease2.Rnw'

###################################################
### code chunk number 1: regsel-heartdisease2.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: regsel-heartdisease2.Rnw:17-18 (eval = FALSE)
###################################################
## source("glmOSCAR_101028.r")


###################################################
### code chunk number 3: regsel-heartdisease2.Rnw:21-43 (eval = FALSE)
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
### code chunk number 4: regsel-heartdisease2.Rnw:49-53 (eval = FALSE)
###################################################
## ######### c fixed
## c.seq<-0.2
## t.seq<-seq(0.01,0.99,length=99)
## oscarR<-glm.oscar(y,X,family,t.seq=t.seq,c=c.seq,epsilon=1e-8)


###################################################
### code chunk number 5: regsel-heartdisease2.Rnw:56-57 (eval = FALSE)
###################################################
## Path<-oscarR$Beta.std[,-1]


###################################################
### code chunk number 6: regsel-heartdisease2.Rnw:60-65 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path*sqrt(n),type="l",
## ylab=ylab.text,xlab=xlab.text,main="oscarR (c=0.2)",lwd=size.width, col=colour)
## axis(4, at = Path[99, ]*sqrt(n), labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


###################################################
### code chunk number 7: regsel-heartdisease2.Rnw:68-71 (eval = FALSE)
###################################################
## c.seq<-0.5
## t.seq<-seq(0.01,0.99,length=99)
## oscarR<-glm.oscar(y,X,family,t.seq=t.seq,c=c.seq,epsilon=1e-8)


###################################################
### code chunk number 8: regsel-heartdisease2.Rnw:74-75 (eval = FALSE)
###################################################
## Path<-oscarR$Beta.std[,-1]


###################################################
### code chunk number 9: regsel-heartdisease2.Rnw:77-82 (eval = FALSE)
###################################################
## par(oma=oma.vec,cex.axis=size.axis,cex.lab=size.axis,cex.main=size.main)
## matplot(rowSums(abs(Path))/max(rowSums(abs(Path))),Path*sqrt(n),type="l",
## ylab=ylab.text,xlab=xlab.text,main="oscarR (c=0.5)",lwd=size.width, col=colour)
## axis(4, at = Path[99, ]*sqrt(n), labels = colnames(X), adj = 0, las = 1,
## cex.axis=size.right)


