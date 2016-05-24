data.weibull <-
function(n,nco){
data.w <- matrix(runif(n*nco),ncol=nco,nrow=n)
colnames(data.w) <- paste("x",1:nco,sep="") 
delta <- rbinom(n,1,0.75)
intercept <- matrix(1,ncol=1,nrow=n)
X <- cbind(intercept,data.w)
Xbeta <- X%*%matrix(1,ncol=1,nrow=ncol(X))
ftime <- NULL
for( j in 1:n){
ftime[j] <- rweibull(1,shape=2,scale=Xbeta[j])
}
wbdat <- data.frame(ftime,data.w,delta)
return(wbdat)
}
