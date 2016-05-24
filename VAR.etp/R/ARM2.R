ARM2 <-
function(x,y,p){
n=length(y); k = ncol(x)
ymat <- matrix(y[(p+1):n])
xmat0 <- matrix(1,nrow=n-p,k*p+1)
index1 = p:(n-1); 
if (p==1) index2=2:(k+1)
if (p>1) index2=seq.int(2,k*p,p)
for (i in 1:p){
xmat0[,index2] <- x[index1,];
index1=index1-1; index2=index2+1
}
b <- solve(crossprod(xmat0)) %*% crossprod(xmat0,ymat)
e <- ymat - xmat0 %*% b
s2 <- sum(e^2)/(nrow(xmat0)-length(b))
covb <- s2*solve(crossprod(xmat0))

# Estimation of k predctors, AR(p) model
amat = matrix(NA,nrow=p+1,ncol=k)
covamat = matrix(NA,nrow=p,ncol=p*k)
emat=matrix(NA,nrow=n-p,k)
x_mat = matrix(NA,nrow=p+1,ncol=k)
index = 1:p; 
for (i in 1:k){
A <- OLS.AR(x[,i],p)
amat[,i] <- A$coef; covamat[,index] <- A$covmat[1:p,1:p]; emat[,i] <- A$resid
tem1=A$xmat; tem2=A$resid
x_mat[,i]=solve(t(tem1) %*% tem1) %*% t(tem1) %*% tem2
index=index+p; 
}

# Calculation of AIC and BIC
mat1 <- cbind(e,emat); mat <- (t(mat1) %*% mat1)/(n-p)
aic = log(det(mat)) + 2* ((k+1)*p)/(n-p)
bic = log(det(mat)) + log(n-p)* ((k+1)*p)/(n-p)

return(list(aic=aic,bic=bic))
}
