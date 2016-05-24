PPrect <-
function(n,l) {
M<-3
s<-n*l;a<-c();b<-c();c<-c()
A<-matrix(1:s, ncol = l, byrow=TRUE)
for (i in 1:n) {
for (j in 1:l) {
a1<-A[i,];b2<-A[-i,];b1<-b2[,j]
c1<-A[-i,];c1<-c1[,-j];c2<-as.vector(c1);f<-length(c2)
a<-c(a,a1);b<-c(b,b1);c<-c(c,c2)}}
AA<-matrix(a, ncol = l, byrow=TRUE)
BB<-matrix(b, ncol = n-1, byrow=TRUE)
CC<-matrix(c ,ncol = f, byrow=TRUE)
Q<-Reduce("cbind",list(AA,"||",BB,"||",CC))
k1<-l
k2<-(n-1)
k3<-(n-1)*(l-1)
lam1<-l+(n-1)*(l-2)
lam2<-(n-2)+(n-2)*(l-1)
lam3<-(n-2)*(l-2)
return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,K=c(k1,k2,k3),lamda=c(lam1,lam2,lam3)))}
