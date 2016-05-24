PPdiv <-
function(n,l){
M<-2
BB<-NULL;AA<-NULL;s<-n*l;A<-matrix(1:s, ncol = l, byrow=TRUE)
for (i in 1:n) {
G1<-t(matrix(rep(A[i,],l),ncol=l))
B<-A[-i,]; B<-as.vector(B)
B1<-t(matrix(rep(B,l),ncol=l))
BB<-rbind(BB,B1);AA<-rbind(AA,G1)}
Q<-Reduce("cbind",list(AA,"||",BB))
k1<-l
k2<-(n-1)*l
lam1<-l+l*(n-1)
lam2<-(n-2)*l
return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,K=c(k1,k2),lamda=c(lam1,lam2)))}
