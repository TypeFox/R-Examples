BIB <-
function(m){
p=2;lin<-(p^(m+1));coll<-m+1;prg<-matrix(nrow=lin, ncol=coll)
for (j in 1:coll ) {prg[,j]<-c(rep(0,lin/(p^j)),rep(1,lin/(p^j)))}
prg<-prg[-1,];t<-dim(prg)[1];f<-dim(prg)[2]
fill<-function(i,mat) {x<-c();m<-c();w<-which(mat[i,]==1)
for (j in 1:t) {
for (k in w) {x<-c(x,prg[j,k])}
if (sum(x)%%2==0) {m<-c(m,j)}
x<-c()}
m}
k<-(t-1)/2;bie<-matrix(nrow=t, ncol=k)
for (i in 1:t) {bie[i,]<-fill(i,prg)}
v <- sort(unique(as.vector(bie)))
V<-length(v);T<-bie[1,1];R<-length(which(T==bie))
B<-dim(bie)[1];K<-dim(bie)[2]
return(list(V=V,B=B,R=R,K=K,BIB=bie))}
