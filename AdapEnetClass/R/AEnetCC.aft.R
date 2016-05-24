AEnetCC.aft <-
function(X, Y, delta, weight, C=1,s = 1, lambda2)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")
weight[weight==0]<-0.001
w <-1/weight
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)

Xextra<-diag(sqrt(lambda2), p)
Yextra<-c(rep(0,p))
delta.extra<-c(rep(1,p))

X<-rbind(X, Xextra)
for (j in 1:p)
X[,j] <-X[,j]/w[j]
Y<-c(Y, Yextra)
delta<-c(delta, delta.extra)
meanx <- apply(X,2,mean)
normx <- sqrt(apply(X^2,2,sum))
meany <- mean(Y)
Y <- Y-meany
X <- t(t(X)/normx)

ncensor <- n+p - sum(delta)
nB <- 2 * p + ncensor
Dmat <- diag(C, nB)
DTmp <- cbind(X[delta==1, ], -X[delta==1, ])
Dmat[1:(2*p), 1:(2*p)] <- t(DTmp) %*% DTmp + diag(lambda2, 2 * p)
dvec <- c( drop(Y[delta==1]) %*% DTmp, numeric(ncensor) )
Aind <- matrix(0, 2*p + 2, nB+1)
Aind[1, 1:ncensor] <- rep(2*p+1, ncensor)
Aind[1, (ncensor+1):nB] <- rep(1, 2*p)
Aind[1, (nB+1)] <- 2*p
for(i in 1:(2*p))
Aind[(i+1), 1:ncensor] <- rep(i, ncensor)
Aind[(2*p+2), 1:ncensor] <- seq(2*p + 1, nB)
Aind[2, (ncensor+1):(ncensor+2*p)] <- seq(1, 2*p)
Aind[2:(2*p+1), nB+1] <- seq(1, 2*p)
Amat <- cbind( rbind(t(X[delta==0,]), t(-X[delta==0,]),
rep(1, ncensor)),
rbind(rep(1, 2*p), matrix(0, 2*p, 2*p)),
c(rep(-1, 2*p), 0) )
bvec <- c(Y[delta==0], numeric(2*p), -s)
qpobj <- solve.QP.compact(Dmat,dvec,Amat,Aind,bvec)
beta <- (qpobj$solution)[1:p] - (qpobj$solution)[(p+1):(2*p)]
beta=t(t(beta)/normx)
beta <-(1+lambda2)*(beta/w)
return( c(beta) )
}
