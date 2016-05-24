### R code from vignette source 'SCMRvin.Rnw'

###################################################
### code chunk number 1: SCMRvin.Rnw:40-52
###################################################
library("SimCorMultRes")
set.seed(1)
N <- 500
ncategories <- 4
clsize <- 3
Xmat <- matrix(rnorm(N),N,ncategories)
betas <- c(1,2,3,4,5,6)
lin.pred <- matrix(c(betas[c(2,4,6)],0),N,4,byrow=TRUE)*Xmat+
           matrix(c(betas[c(1,3,5)],0),N,4,byrow=TRUE)
lin.pred <- matrix(lin.pred,N,ncategories*clsize)
cor.matrix <- diag(1,12)
Y <- rmult.bcl(clsize,ncategories,lin.pred,cor.matrix)


###################################################
### code chunk number 2: SCMRvin.Rnw:55-56
###################################################
head(Y$Ysim)


###################################################
### code chunk number 3: SCMRvin.Rnw:87-94
###################################################
set.seed(12345)
N <- 500
clsize <- 4
intercepts <- c(-1.5,-0.5,0.5,1.5)
cor.matrix <- toeplitz(c(1,0.85,0.5,0.15))
lin.pred <- rsmvnorm(N,toeplitz(c(1,rep(0.85,clsize-1))))
Y <- rmult.clm(clsize,lin.pred,cor.matrix,intercepts,"probit")


###################################################
### code chunk number 4: SCMRvin.Rnw:97-98
###################################################
head(Y$Ysim)


###################################################
### code chunk number 5: SCMRvin.Rnw:123-131
###################################################
set.seed(1)
N <- 500
clsize <- 4
intercepts <- c(-1.5,-0.5,0.5,1.5)
cor.matrix <- diag(1,16)
x <- rnorm(N)
lin.pred <- matrix(rep(x,clsize),N,clsize,byrow=TRUE)
Y <- rmult.crm(clsize,lin.pred,cor.matrix,intercepts,link="probit")


###################################################
### code chunk number 6: SCMRvin.Rnw:134-135
###################################################
head(Y$Ysim)


###################################################
### code chunk number 7: SCMRvin.Rnw:162-169
###################################################
set.seed(1)
N <- 500
clsize <- 4
intercepts <- 1
cor.matrix <- toeplitz(c(1,0.85,0.5,0.15))
lin.pred <- matrix(rnorm(N),N,clsize)
Y <- rbin(clsize,lin.pred,cor.matrix,intercepts,"probit")   


###################################################
### code chunk number 8: SCMRvin.Rnw:172-173
###################################################
head(Y$Ysim)


