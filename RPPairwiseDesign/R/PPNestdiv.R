PPNestdiv <-
function(n,l,w){
if (w <= 2 ) {"Choose w > 2"}
else
M<-3
s<-l*n;A<-NULL
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A)
N<-n*w; L<-l*l
G1<-matrix(nrow=N, ncol=L)
G2<-matrix(nrow=N, ncol=L*(n-1))
for (i in 1:N) {
G1[i,]<-rep(B[i,],l)
I<-((i-1)%/%n)+1
J<-(i%%n)
if (J==0){J<-n}
X<-as.vector(t(A[[I]][-J,]))
G2[i,]<-rep(X,l)}
NN<-n*w*l; LL<-l*n*(w-1)
G3<-NULL
for (o in 1:w) {
Y<-Reduce("rbind",A[-o])
Ys<-as.vector(t(Y));tt<-rep(t(Ys),n*l)
t3<-matrix(tt, ncol=n*l*(w-1),byrow=TRUE)
G3[[o]]<-t3}
T1<-matrix(t(G1), ncol=l,byrow=TRUE)
T2<-matrix(t(G2), ncol=l*(n-1),byrow=TRUE)
T3<-Reduce("rbind",G3)

Q<-Reduce("cbind",list(T1,"||",T2,"||",T3))

k1<-l
k2<-l*(n-1)
k3<-(w-1)*l*n
lam1<-l+l*(n-1)+(w-1)*l*n
lam2<-(n-2)*l+(w-1)*l*n
lam3<-(w-2)*l*n

return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,k=c(k1,k2,k3),lamda=c(lam1,lam2,lam3)))}
