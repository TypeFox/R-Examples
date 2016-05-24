TEST2 <-
function(f){
TEST <-
function(Nt){
QB<-0
QS<-0
QR<-0
PB<-0
PR<-0
PS<-0
ZB<-0
ZS<-0
ZR<-0
V<-matrix(0,3,3)
d<-matrix(0,3,1)
for(i in 1:4){
for(j in 1:4){
if((i<j))
if(Nt[i,j]+Nt[j,i]==0){
}
else{
QB<-QB+(((Nt[i,j]-Nt[j,i])^2)/(Nt[i,j]+Nt[j,i]))}
}
}
for(i in 1:3){
d[i]<-(sum(Nt[i,])-sum(Nt[,i]))
for(j in 1:3)
if(i==j)V[i,j]<-(sum(Nt[i,])+sum(Nt[,i])-2*sum(Nt[i,i]))
else V[i,j]<--1*(Nt[i,j]+Nt[j,i])
}
QS<-t(d)%*%(solve(V))%*%d
QR<-QB-QS
PB<-(1-pchisq(QB,6))
Z<-qnorm(PB)
PS<-(1-pchisq(QS,3))
ZS<-qnorm(PS)
PR<-(1-pchisq(QR,3))
ZR<-qnorm(PR)
r<-matrix(c(PB,PS,PR),3,1)
r
}
n<-length(dim(f))
B1<-NULL
S1<-NULL
I1<-NULL
for(i in 1:n){
B<-matrix(0,n,1)
S<-matrix(0,n,1)
I<-matrix(0,n,1)
for(j in 1:n){
if(j!=i){
k<-apply(f,c(i,j),sum)
testr<-TEST(k)
}
else{
testr<-matrix(0,3,1)
}
B[j]<-testr[1]
S[j]<-testr[2]
I[j]<-testr[3]
}
B1<-cbind(B1,B)
S1<-cbind(S1,S)
I1<-cbind(I1,I)
}
list(Bowker.test=as.dist(round(B1,4)),
Stuart.test=as.dist(round(S1,4)),
Internal.test=as.dist(round(I1,4)))
}
