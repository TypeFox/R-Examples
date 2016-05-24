mvtmeta_re <-
function(y, cov) {
fixed<-mvtmeta_fe(y, cov)
p<-dim(y)[1]
n<-dim(y)[2]
yres<-y-fixed$beta
cov_i<-array(rep(NA, p*p*n), c(p, p, n))
cov_i_yressq<-array(rep(NA, p*p*n), c(p, p, n))
tmpphi<-array(rep(NA, p*p*n), c(p, p, n))
tmpy<-matrix(rep(NA, p*n), p, n)
for(i in 1:n) {
cov_i[,,i]<-solve(cov[,,i])
cov_i_yressq[,,i]<-cov_i[,,i]%*%yres[,i]%*%t(yres[,i])
tmpphi[,,i]<-cov_i[,,i]-cov_i[,,i]%*%fixed$cov%*%cov_i[,,i]
}
phi<-apply(tmpphi, 1:2, sum)
cov_i_yressqsum<-apply(cov_i_yressq, 1:2, sum)
A<-solve(phi,cov_i_yressqsum-diag(rep(n-1, p)))
TT<-(A+t(A))/2
eig<-eigen(TT)
negeigen<-sum(eig$values<0)
if(negeigen==0) Tpd<-TT
else Tpd<-eig$vectors%*%diag(pmax(eig$values, 0))%*%t(eig$vectors)
for(i in 1:n) {
cov[,,i]<-cov[,,i]+Tpd
cov_i[,,i]<-solve(cov[,,i])
tmpy[,i]<-cov_i[,,i]%*%y[,i]
}
covbeta_i<-apply(cov_i, 1:2, sum)
covbeta<-solve(covbeta_i)
beta<-as.vector(covbeta%*%apply(tmpy, 1, sum))
return(list(beta=beta, cov=covbeta, between=Tpd, negeigen=negeigen))
}

