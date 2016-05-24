pem<-function(x){
n <- sum(x)
R<-dim(x)[1]
C<-dim(x)[2]
sr <- rowSums(x)
SR<-matrix(rep(sr,C),nrow=R)
sc <- colSums(x)
SC<-matrix(rep(sc,R),nrow=R,byrow=TRUE)
E <- outer(sr, sc, "*")/n
denom<-pmin(SR,SC)-E
Min<-pmax(SR+SC-n,0)
denom[x-E<0]<--(Min[x-E<0]-E[x-E<0])
PEM<-100*(x-E)/denom
PEM
}
