Enet.wls <-
function(X, Y, delta)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and delta don't match!")
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)
fit<-glmnet(X, Y)
beta<-fit$beta[,dim(fit$beta)[2]]
return(list(beta=beta, fit=fit))
}
