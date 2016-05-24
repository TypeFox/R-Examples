ARM <-
function(x,y,p)
{
n=length(y)
ymat <- matrix(y[(p+1):n])
xmat0 <- matrix(1,nrow=n-p,p+1)
index = p:(n-1)
for (i in 1:p){
xmat0[,i+1] <- x[index];
index=index-1
}
b <- solve(crossprod(xmat0)) %*% crossprod(xmat0,ymat)
e <- ymat - xmat0 %*% b
s2 <- sum(e^2)/(nrow(xmat0)-length(b))
covb <- s2*solve(crossprod(xmat0))
A <- OLS.AR(x,p)
a <- A$coef; cova <- A$covmat[1:p,1:p]; e1 <- A$resid
 
mat1 <- cbind(e,e1); mat <- (t(mat1) %*% mat1)/(n-p)
aic = log(det(mat)) + 2* (2*p)/(n-p)
bic = log(det(mat)) + log(n-p)* (2*p)/(n-p)
#aic = log( sum(e1^2)/(n-p)) + 2* (p+1)/(n-p)
#bic = log( sum(e1^2)/(n-p)) + log(n)* (p+1)/(n-p)
return(list(b=b,a=a,covb=covb,cova=cova,aic=aic,bic=bic,xmat=xmat0,ymat=ymat))
}
