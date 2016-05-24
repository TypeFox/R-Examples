VAR.Pope <-
function(x,p,type="const")
{
n <- nrow(x); k <- ncol(x)
var1 <- VAR.est(x,p,type)
b <- var1$coef
sigu <- var1$sigu
zz <- var1$zzmat

zmat <- zz[1:(p*k),1:(p*k)]
if( p == 1) 
{gmat <- sigu
amat <- b[,1:(k*p)]}
if( p > 1) 
{gmat <- cbind(rbind(sigu,matrix(0,nrow=k*(p-1),ncol=k)),matrix(0,nrow=k*p,ncol=k*(p-1)))
amat <- rbind(b[,1:(k*p)], cbind(diag(k*(p-1)),matrix(0,nrow=k*(p-1),ncol=k)) )}

lamda <- matrix(eigen(amat)$values)
sum3 <- 0
for(i in 1:length(lamda))
sum3 <- sum3+ lamda[i]*solve(diag(k*p)-lamda[i] * t(amat) )

sum1 <- solve(diag(k*p) - t(amat))
sum2 <- t(amat) %*% solve(diag(k*p)- t(amat) %*% t(amat))
sum4 <- sum1+sum2+sum3

bias <- -(gmat %*% sum4 %*% solve(zmat))/n
bias <- matrix(bias[1:k,],nrow=k); rownames(bias) <- rownames(b); colnames(bias)=colnames(b)[1:(p*k)]

bc <- b[,1:(k*p)]-bias; 
bs <- VAR.adjustP(b,bias,p,type); rownames(bs) <- rownames(b); colnames(bs) <- VAR.names(x,p,type)
es <- VAR.resid(x,bs,var1$zmat,p); colnames(es) <- rownames(b)
sigu <- t(es) %*% es / ( (n-p) -k*p -1)

return(list(coef=Re(bs),resid=Re(es),sigu=Re(sigu),Bias=Re(bias)))
}
