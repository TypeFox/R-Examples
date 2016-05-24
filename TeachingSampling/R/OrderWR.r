OrderWR<-function(N,m,ID=FALSE){
b<-c(1:N)
grilla<-function(a){
A<-seq(1:length(a))
unoA <-rep(1,length(A))
B<-seq(1:length(a))
unoB <-rep(1,length(B))
P1<-kronecker(A,unoB)
P2<-kronecker(unoA,B)
grid<-matrix(cbind(P1,P2),ncol=2)
return(grid)
}

if(m==1){
sam<-as.matrix(b)
}

if(m==2){
sam<-grilla(b)
}

if(m>2){
sam<-grilla(b)
for(l in 3:m){
Sam1<-rep(0,l)
for(j in 1:dim(sam)[1]){
for(k in 1:length(b)){
Sam1<-rbind(Sam1,c(sam[j,],b[k]))
} }
sam<-Sam1[-1,]
}
}

if(ID==FALSE){return(sam)}
else{
a<-dim(sam)
val<-matrix(NA,a[1],a[2])
for(ii in 1:(dim(val)[1])){
for(jj in 1:(dim(val)[2])){
val[ii,jj]<-ID[sam[ii,jj]]
}
}
return(val)
}
}