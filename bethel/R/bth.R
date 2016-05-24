`bth` <-
function(S,T,eps=0.1e-09){
S<-as.data.frame(S)
T<-as.data.frame(T)

#error
if (ncol(S)<=5) stop("ncol(S) must be >5!See the online help: help('bth').")
if (ncol(S)-5!=nrow(T)) stop("ncol(S)-5 must be equal to nrow(S)! See the online help: help('bth').")
if (sum(S[,2:(ncol(S)-4)])==0) stop("All the strata have null variance. Then the optimal sample size is n=1 in each strata!")
if ((t(S[,ncol(S)-3]>0))%*%matrix(rep(TRUE,nrow(S)))!=nrow(S)) stop("The number of units in each strata can not be <=0!")
if ((t(S[,ncol(S)-2]>0))%*%matrix(rep(TRUE,nrow(S)))!=nrow(S)) stop("Sample costs in each strata can not be <=0!")
if ((t(S[,ncol(S)-1]>0))%*%matrix(rep(TRUE,nrow(S)))!=nrow(S)) stop("The minimum sample size can not be <=0!")
if ((t(S[,ncol(S)]>0))%*%matrix(rep(TRUE,nrow(S)))!=nrow(S)) stop("The minimum % sample size can not be <=0!")

#datasets
num<-t(S[,2:(ncol(S)-4)]*(S[,(ncol(S)-3)]^2))
den<-T[,1]^2*T[,2]^2+t(as.matrix(S[,2:(ncol(S)-4)]))%*%as.matrix(S[,(ncol(S)-3)])
if (ncol(S)<=6){ 
A<-matrix(num / as.vector(den))
}
else {A<-t(diag(as.vector(1/den))%*%num)}
c<-cbind(S[,ncol(S)-2])

#bethel
alfa_1<-matrix(rep(1/ncol(A),ncol(A)))
ck_alfa<-matrix(rep(TRUE,ncol(A)))
while (t(ck_alfa)%*%ck_alfa>0){
alfa_0<-alfa_1
x<-sqrt(c)/( sqrt(A%*%alfa_0) %*% sqrt(t(c))%*%sqrt(A%*%alfa_0))
for (i in 1:nrow(x))if (is.infinite(x[i,1])) x[i,1]<-1e+09
alfa_1<-alfa_0 * (t(A)%*%x)^2 / sum(alfa_0 * (t(A)%*%x)^2)
ck_alfa<-abs(alfa_1-alfa_0)>eps
}

#output 
B<-data.frame(S[,1],ceiling(1/x),pmin(pmax(S[,(ncol(S)-1)],ceiling(S[,ncol(S)]*S[,(ncol(S)-3)]),ceiling(1/x)),S[,ncol(S)-3]))
names(B)=c("strata","numBethel","numBethel2")
B
}