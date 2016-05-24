TEST3 <-
function(Farray){
s1<-sum(Farray)
Farray<-Farray/s1
v<-NULL
m<-length(dim(Farray))
r<-dim(Farray)[1]
n<-NULL
one<-rep(1,r*(m-1))
J<-matrix(0,4*(m-1),4*(m-1))
for(i in 1:(m-1)){
J[(r*(i-1)+1):(r*(i-1)+r),(r*(i-1)+1):(r*(i-1)+r)]<-1
}
L1<-t(matrix(t(rep(diag(c(1,1,1,1)),(m-1))),r,(r*(m-1))))
L<-cbind(L1,-(diag(one)))
for(i in 1:m){
vc1<-NULL
for(k in 1:m){
vc<-NULL
if(i==k){
fi<-apply(Farray,i,sum)
vc<-diag(fi)
n1<-matrix(fi,r,1)
n<-rbind(n,n1)
}
if(i<k){
fij<-apply(Farray,c(i,k),sum)
fr<-apply(fij,1,sum)
fc<-apply(fij,2,sum)
vc<-fij
}
if(i>k){
fij<-apply(Farray,c(i,k),sum)
fr<-apply(fij,1,sum)
fc<-apply(fij,2,sum)
vc<-fij
}
vc1<-cbind(vc1,vc)
}
v<-rbind(v,vc1)
}
Ts<-(t(L%*%n)%*%solve((L%*%v%*%t(L))+J)%*%(L%*%n))*s1
Ts
}
