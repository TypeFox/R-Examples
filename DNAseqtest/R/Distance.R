Distance <-
function(F4){
L1<-length(dim(F4))
Dis<-matrix(0,L1,L1)
for(i in 1:L1){
for(j in 1:L1){
if(i!=j){
fij<-apply(F4,c(i,j),sum)
fi<-apply(fij,1,sum)
fj<-apply(fij,2,sum)
Dis[i,j]<--log(det(diag((fi)^-0.5)%*%fij%*%diag((fj)^-0.5)))
}
else{
Dis[i,j]<-0
}
}
}
Dis
}
