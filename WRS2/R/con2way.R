con2way<-function(J,K){
#
# For a  J by K ANOVA design
# create the contrast coefficients for
# doing all pairwise comparisons of
# main effects for Factor A and B and all interactions
#
JK <- J * K
Ja<-(J^2-J)/2
Ka<-(K^2-K)/2
JK<-J*K
conA<-matrix(0,nrow=JK,ncol=Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,]<-1
mat[jj,]<-0-1
conA[,ic]<-t(mat)
}}}
conB<-matrix(0,nrow=JK,ncol=Ka)
ic<-0
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[,k]<-1
mat[,kk]<-0-1
conB[,ic]<-t(mat)
}}}
conAB<-matrix(0,nrow=JK,ncol=Ka*Ja)
ic<-0
for(j in 1:J){
for(jj in 1:J){
if(j < jj){
for(k in 1:K){
for(kk in 1:K){
if(k<kk){
ic<-ic+1
mat<-matrix(0,nrow=J,ncol=K)
mat[j,k]<-1
mat[j,kk]<-0-1
mat[jj,k]<-0-1
mat[jj,kk]<-1
}
conAB[,ic]<-t(mat)
}}}}}
list(conA=conA,conB=conB,conAB=conAB)
}
