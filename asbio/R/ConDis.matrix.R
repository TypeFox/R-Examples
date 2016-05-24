ConDis.matrix<-function(Y1,Y2){
n<-length(Y1)
m<-as.data.frame(matrix(nrow=n,ncol=n,dimnames=list(seq(1,n), seq(1,n))))
for(i in 1:n){
for(j in 1:n){
m[j,i]<-ifelse((Y1[i]>Y1[j]&Y2[i]>Y2[j])|(Y1[i]<Y1[j]&Y2[i]<Y2[j]),1,
ifelse((Y1[i]<Y1[j]&Y2[i]>Y2[j])|(Y1[i]>Y1[j]&Y2[i]<Y2[j]),-1,0))
}}
m[upper.tri(m,diag=TRUE)]<-NA
m
}
