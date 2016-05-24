data.weibull.reg <-
function(n,ncop,shape){
mu <- rep(5,ncop)
Sigma <- matrix(0.7,nrow=ncop,ncol=ncop)+diag(ncop)*0.25
data.wr <- mvrnorm(n, mu,Sigma)
colnames(data.wr) <- paste("x",1:ncop,sep="")
delta <- rbinom(n,1,0.75)
intercept <- matrix(1,ncol=1,nrow=n)
X <- cbind(intercept,data.wr)
Xbeta <- as.matrix(X)%*%matrix(1,ncol=1,nrow=ncol(X))
ftime <- NULL
for( j in 1:n){
ftime[j] <- rweibull(1,shape=shape,scale=Xbeta[j])
}
wbdat <- data.frame(ftime,data.wr,delta)
return(wbdat)
}
