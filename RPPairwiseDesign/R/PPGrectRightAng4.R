PPGrectRightAng4 <-
function(n,l,w) {
m<-w
if (l < 3) {"Choose l > 2"}
else {
M<-4
s<-l*n;A<-NULL
for (i in 1:m){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s
A[[i]]<-A[[i]]+z}
B<-Reduce("rbind",A)
b<-Reduce("cbind",A)
N<-n*m; L<-l*l
G3<-NULL;G1<-matrix(nrow=N, ncol=L)
G2<-NULL
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
NN<-n*m*l; LL<-l*n*(m-1)
G4<-NULL
for (o in 1:m) {
Y<-Reduce("rbind",A[-o])
Ys<-as.vector(t(Y))
tt<-rep(t(Ys),n*l)
t4<-matrix(tt, ncol=n*l*(m-1),byrow=TRUE)
G4[[o]]<-t4}
d<-dim(b)[1]
for (p in 1:d) {
b2<-b[-p,];G<-t(b2);G2<-rbind(G2,G)}
T1<-matrix(t(G1), ncol=l,byrow=TRUE);T2<-G2;T3<-G3
T4<-Reduce("rbind",G4)
Q<-Reduce("cbind",list(T1,"||",T2,"||",T3,"||",T4))
k1<-dim(T1)[2]
k2<-dim(T2)[2]
k3<-dim(T3)[2]
k4<-dim(T4)[2]
lam1<-l+(l-2)*(n-1)+(w-1)*l*n
lam2<-(n-2)+(l-1)*(n-2)+(w-1)*l*n
lam3<-(l-2)*(n-2)+(w-1)*l*n
lam4<-(w-2)*l*n
return(list(RPPBD=noquote(Q),v=s,b=s*M,r=s,k=c(k1,k2,k3,k4),lamda=c(lam1,lam2,lam3,lam4)))}}
