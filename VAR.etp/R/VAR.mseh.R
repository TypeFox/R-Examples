VAR.mseh <-
function(x,p,h,type) {
n <- nrow(x); k <- ncol(x)
var1 <- VAR.est(x,p,type)
b <- var1$coef[,1:(k*p+1)]
zz <- var1$zzmat[1:(k*p+1),1:(k*p+1)]
sigu <- var1$sigu

mf <- VAR.mainf(b,p,h)

if (h == 1) sigyh <- ( ((n-p)+k*p+1)/(n-p) ) * sigu

if (h > 1) {
if (p == 1) bb <- rbind(b,cbind(matrix(0,ncol=k*p),1))
if (p > 1)  bb <- rbind(b, cbind(diag(k*(p-1)),matrix(0,nrow=k*(p-1),ncol=k),matrix(0,nrow=k*(p-1),ncol=1)),cbind(matrix(0,ncol=k*p),1))

b0 <- diag(k*p+1); bbq <- b0
for( i in 1:(h-1)){
b0 <- b0 %*% bb
bbq <- cbind(bbq,b0)}

index1 <- seq((h-1)*(k*p+1)+1,length.out=p*k+1); index3 <- 1:k
sum1 <- matrix(0,nrow=k,ncol=k); sigy <- sum1

for (i in 0:(h-1))
   {
   sigy <- sigy + mf[,index3]%*%sigu%*% t(mf[,index3])
   index2 <- seq((h-1)*(k*p+1)+1,length.out=p*k+1); index4 <- 1:k
   for (j in 0:(h-1)){
   tem <- sum( diag( t(bbq[,index1]) %*% solve(zz) %*% bbq[,index2] %*% zz) ) * mf[,index3] %*% sigu %*% t(mf[,index4])
   sum1 <- sum1 + tem
   index2 <- index2 - (p*k+1); index4 <- index4 + k;}
   index1 <- index1-(p*k+1); index3 <- index3+k
   }
sigyh <- sigy + sum1/(n-p);
}
return(sigyh)
}
