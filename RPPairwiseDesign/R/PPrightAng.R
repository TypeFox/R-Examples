PPrightAng <-
function(n,l,w) {
m<-w
if (l < 3 ) {"Choose l > 2"}
else {
M<-4
s<-l*n;A<-NULL
for (i in 1:m){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s;A[[i]]<-A[[i]]+z}
B<-Reduce("rbind",A)
N<-n*m; L<-l*l
G1<-matrix(nrow=N, ncol=L);G3<-matrix(nrow=N, ncol=l*(m-1)*l)
G2<-matrix(nrow=N, ncol=l*l*(n-1));G4<-matrix(nrow=N, ncol=l*l*(n-1)*(m-1))
for (i in 1:N) {
G1[i,]<-rep(B[i,],l)
I<-((i-1)%/%n)+1
J<-(i%%n)
if (J==0){J<-n}
X<-as.vector(t(A[[I]][-J,]))
G2[i,]<-rep(X,l)
AA<-A[-I]
H<-Reduce("cbind",AA)
Y<-H[J,]
Z<-as.vector(t(H[-J,]))
G3[i,]<-rep(Y,l);G4[i,]<-rep(Z,l)}
T1<-matrix(t(G1), ncol=l,byrow=TRUE)
T3<-matrix(t(G3), ncol=l*(m-1),byrow=TRUE)
T2<-matrix(t(G2), ncol=l*(n-1),byrow=TRUE)
T4<-matrix(t(G4), ncol=l*(n-1)*(m-1),byrow=TRUE)
Q<-Reduce("cbind",list(T1,"||",T2,"||",T3,"||",T4))
k1<-dim(T1)[2]
k2<-dim(T2)[2]
k3<-dim(T3)[2]
k4<-dim(T4)[2]
lam1<-l+(n-1)*l+(w-1)*l+(w-1)*(n-1)*l
lam2<-(n-2)*l+(w-1)*l*(n-2)
lam3<-(w-2)*l+(w-2)*l*(n-1)
lam4<-(w-2)*(n-2)
return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,K=c(k1,k2,k3,k4),lamda=c(lam1,lam2,lam3,lam4)))}}
