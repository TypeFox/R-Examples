VARB.Pope <-
function(x,p,type="const")
{
n <- nrow(x); k <- ncol(x)
var1 <- VARB.est(x,p,type)
b <- var1$coef
sigu <- var1$sigu
zz <- var1$zzmat

zmat <- zz[1:(p*k),1:(p*k)]
if( p == 1) 
{gmat <- sigu
amat <- b[,1:(k*p)]}
if( p > 1) 
{gmat <- cbind( matrix(0,nrow=k*p,ncol=k*(p-1)),rbind(matrix(0,nrow=k*(p-1),ncol=k), sigu) )
amat <- {}
index <- 1:k
for (i in 1:p){
amat <- cbind(b[,index],amat)
index <- index+k}
amat <- rbind(cbind(matrix(0,nrow=k*(p-1),ncol=k),diag(k*(p-1))),amat )
}

lamda <- matrix(eigen(amat)$values)
sum3 <- 0
for(i in 1:length(lamda))
sum3 <- sum3+ lamda[i]*solve(diag(k*p)-lamda[i] * t(amat) )

sum1 <- solve(diag(k*p) - t(amat))
sum2 <- t(amat) %*% solve(diag(k*p)- t(amat) %*% t(amat))
sum4 <- sum1+sum2+sum3

bias <- -(gmat %*% sum4 %*% solve(zmat))/n
bias <- matrix(bias[1:k,],nrow=k)

bc <- b[,1:(k*p)]-bias; 
bs <- VAR.adjustP(b,bias,p,type); rownames(bs) <- rownames(b); colnames(bs) <- VAR.names(x,p,type)
es <- VAR.resid(x,bs,var1$zmat,p); colnames(es) <- rownames(b)
sigu <- t(es) %*% es / ( (n-p) -k*p -1)

return(list(coef=Re(bs),resid=Re(es),sigu=Re(sigu)))
}
