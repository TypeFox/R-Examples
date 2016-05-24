#source("charge.R")
library("LINselect")

set.seed(15) # pour pouvoir comparer des exécutions sucessives
ex <- simulData(p=25,n=25,r=0.5,rSN=10)

resVARselect <- VARselect(ex$Y,ex$X,exhaustive.dmax=4)
resVARselect$summary

resVARselect <-
  VARselect(ex$Y,ex$X,normalize=FALSE,dmax=15,exhaustive.dmax=4,verbose=FALSE) 
resVARselect$summary

## Time consuming
# X <- ex$X
# X[,1] <- X[,1]/10
# X[,2] <- X[,3]*10
# resVARselect <- VARselect(ex$Y,X,dmax=20,exhaustive.dmax=2,verbose=FALSE)
# resVARselect$summary

# resVARselect <- VARselect(ex$Y,X,normalize=FALSE,dmax=20,exhaustive.dmax=2,verbose=FALSE)
# resVARselect$summary
# 
# resVARselect <- VARselect(ex$Y,ex$X,method="exhaustive",exhaustive.maxdim=5000)
# resVARselect$summary
# 
# resVARselect <- VARselect(ex$Y,ex$X,method="en",en.lambda=5,dmax=20)
# resVARselect$summary
# 
# set.seed(15)
# ex <- simulData(p=100,n=100,r=0.5,rSN=10)
# 
# methodF=c("lasso","ridge","pls","en","ALridge","ALpls","rF")
# resVARselect <- VARselect(ex$Y,ex$X,method=methodF)
# resVARselect$summary
# 
# resVARselect <- VARselect(ex$Y,ex$X,method="rF",rF.lmtry=c(1,2,5))
# resVARselect

# 
# set.seed(15) 
# ex <- simulData(p=100,n=50,r=0.95,rSN=10)
# 
# resVARselect <- VARselect(ex$Y,ex$X,method=methodF)
# resVARselect$summary
# 
# set.seed(15)
# ex <- simulData(p=1000,n=200,r=0.95,rSN=10)
# methodF=c("lasso","ridge","pls","en","ALridge","ALpls")
# resVARselect <- VARselect(ex$Y,ex$X,method=methodF,dmax=50,
#                       en.lambda=c(0.01,1,5),
#                       ridge.lambda=c(0.01,1,5),
#                       ALridge.lambda=c(0.01,1,5))
# resVARselect$summary
# 
# set.seed(15) 
# p <- 5000 
# n <- 100
# X <-  matrix(rnorm(n*p),ncol=p,nrow=n)
# Y <- rnorm(n)
# resVARselect <- VARselect(Y,X,method="ALpls",ALpls.ncomp=4,dmax=50)
# resVARselect







