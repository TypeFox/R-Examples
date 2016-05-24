Gen <-
function(n,mat){
A<-mat[n,];B<-mat[-n,]
a<-dim(B)[1];b<-dim(B)[2]
e<-(b-1)/2
x1<-matrix(nrow=a, ncol=e)
for (j in 1:a) {for (i in 1:b) {if (all (B[j,i]!=A)) {B[j,i]<-0}}}
for (i in 1:a) {x1[i,]<-B[i,][B[i,]>0]} 
x1<-unique(x1)
v <- sort(unique(as.vector(x1)))
V<-length(v);T<-x1[1,1];R<-length(which(T==x1))
B<-dim(x1)[1];K<-dim(x1)[2]
return(list(V=V,B=B,R=R,K=K,BIB2=x1))}
