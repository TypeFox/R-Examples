Pt <-
function(S,Pi,t){
ax<-matrix(0,nrow=4,ncol=4)
fx<-matrix(0,nrow=4,ncol=4)
Rx<-S%*%Pi
a1<-eigen((sqrt(Pi))%*%Rx%*%(solve(sqrt(Pi))),symmetric=T)$values
a2<-eigen((sqrt(Pi))%*%Rx%*%(solve(sqrt(Pi))),symmetric=T)$vectors
 for(j in 1:4){
 ax<-exp(a1[j]*(t))*a2[,j]%*%t(a2[,j])
 fx<-fx+ax
 }
ptx<-(solve(sqrt(Pi)))%*%fx%*%(sqrt(Pi))
ptx
}
