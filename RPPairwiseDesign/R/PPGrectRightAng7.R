PPGrectRightAng7 <-
function(n,l,w) {
m<-w
if (l < 3 | m < 3 | n < 3) {"Choose n > 2 & l > 2 & w > 2"}
else {
M<-7
s<-l*n;A<-NULL
for (i in 1:m){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s
A[[i]]<-A[[i]]+z}
B<-Reduce("rbind",A)
b<-Reduce("cbind",A);N2<-n*l
N<-n*m; L<-l*l;G3<-NULL
G1<-matrix(nrow=N, ncol=L);G2<-NULL
G4<-NULL;G5<-NULL;G6<-NULL;G7<-NULL
M<-m-1
for (i in 1:N) {
G1[i,]<-rep(B[i,],l)
I<-((i-1)%/%n)+1
J<-(i%%n)
if (J==0){J<-n}
x<-A[[I]][-J,]
X<-matrix(nrow=l, ncol=(n-1)*(l-1))
for (k in 1:l){
xx<-x[,-k]
X[k,]<-as.vector(t(xx))}
G3<-rbind(G3,X)}
d<-dim(b)[1]
for (p in 1:d) {
b2<-b[-p,]
G<-t(b2)
G2<-rbind(G2,G)}
tt<-l*(m-1)
for (v in 1:m){
BB<-A[-v];B2<-Reduce("cbind",BB);B3<-Reduce("rbind",BB)
BD<-NULL
for (k in 1:l){
ve<-seq(k,tt,l);GG<-B2[,-ve];G5<-rbind(G5,GG)}
for (k in 1:n) {
vs<-seq(k,tt,n);FF<-B3[-vs,]
gG<-t(FF);G6<-rbind(G6,gG)
for (o in 1:l) {
BD<-B3[-vs,-o];be<-as.vector(BD)
G7<-rbind(G7,be)}}}
g4<-matrix(nrow=N2, ncol=m)
for (h in 1:m){
g4[,h]<-as.vector(t(A[[h]]))}
for (H in 1:m){
g<-g4[,-H];G4<-rbind(G4,g)}
T1<-matrix(t(G1), ncol=l,byrow=TRUE);T2<-G2
T3<-G3;T4<-G4;T5<-G5;T6<-G6
rownames(G7)<-numeric(0);T7<-G7
Q<-Reduce("cbind",list(T1,"||",T2,"||",T3,"||",T4,"||",T5,"||",T6,"||",T7))
k1<-dim(T1)[2]
k2<-dim(T2)[2]
k3<-dim(T3)[2]
k4<-dim(T4)[2]
k5<-dim(T5)[2]
k6<-dim(T6)[2]
k7<-dim(T7)[2]
lam1<-l+(l-2)*(n-1)+(w-1)*(l-2)+(w-1)*(n-1)*(l-2)
lam2<-(n-2)+(l-1)*(n-2)+(w-1)*(n-2)+(w-1)*(l-1)*(n-2)
lam3<-(l-2)*(n-2)+(w-1)*(l-2)*(n-2)
lam4<-(w-2)+(w-2)*(l-1)+(w-2)*(n-1)+(w-2)*(l-1)*(n-1)
lam5<-(w-2)*(l-2)+(w-2)*(l-2)*(n-1)
lam6<-(w-2)*(n-2)+(w-2)*(l-1)*(n-2)
lam7<-(w-2)*(n-2)*(l-2)
return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,k=c(k1,k2,k3,k4,k5,k6,k7),lamda=c(lam1,lam2,lam3,lam4,lam5,lam6,lam7)))}}
