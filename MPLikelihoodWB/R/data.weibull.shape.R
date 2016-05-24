data.weibull.shape <-
function(n,nco,shape){
data.ws <- matrix(runif(n*nco),ncol=nco,nrow=n)
colnames(data.ws) <- paste("x",1:nco,sep="") 

delta <- rbinom(n,1,0.75)
intercept <- matrix(1,ncol=1,nrow=n)
X <- cbind(intercept,data.ws)
Xbeta <- X%*%matrix(1,ncol=1,nrow=ncol(X))
ftime <- NULL
for( j in 1:n){
ftime[j] <- rweibull(1,shape=shape,scale=Xbeta[j])
}

wbdat <- data.frame(ftime,data.ws,delta)
return(wbdat)
}
