.packageName <- "sparseLDA"

normalize <- function(X){
# normalize the columns of X to have zero mean and length one.
# output is:
# $Xc : Normalized X
# $mx : Mean of X
# $vx : Length of columns of X
# $Id : Logical indicating which columns have a length different from zero
 
n <- dim(X)[1]
p <- dim(X)[2]
mx <- apply(X,2,mean)
m <- rep(mx,n)
dim(m) <- c(p,n)
X <- X-t(m)
vx <- sqrt(apply(X^2,2,sum))
Id <- vx!=0
v <- rep(vx,n)
dim(v) <- c(p,n)
v <- t(v)
X <- X[,Id]/v[,Id]
list(Xc=X,mx=mx,vx=vx,Id=Id)
}


