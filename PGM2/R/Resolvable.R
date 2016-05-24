Resolvable <-
function(n,mat) {B<-mat;C<-B[n,];B<-B[-n,]
a<-dim(B)[1]
c<-dim(B)[2];e<-(c+1)/2
for (j in 1:a) {for (i in 1:c) {if (any (B[j,i]==C)) {B[j,i]<-0}}}
X<-matrix(nrow=a, ncol=e)
for (i in 1:a) {X[i,]<-B[i,][B[i,]>0]}
v <- sort(unique(as.vector(X)))
V<-length(v);T<-X[1,1];R<-length(which(T==X))
B<-dim(X)[1];K<-dim(X)[2]
return(list(V=V,B=B,R=R,K=K,RBIB=X))}
