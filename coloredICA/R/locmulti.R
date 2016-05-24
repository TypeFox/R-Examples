locmulti <-
function(x0,l_period,n,freq,h){
X=cbind(rep(1,n),freq[,1]-x0[1],freq[,2]-x0[2])
K=diag(kern(x0,h,freq)$v)
a=solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%l_period
list(ahat=a)
}
