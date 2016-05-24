kron<-function(m1,m2){
#  compute the Kronecker product of the two matrices m1 and m2.
#
m1<-as.matrix(m1) # Vectors of length p are converted to a p by 1 matrix
m2<-as.matrix(m2)
kron<-vector(mode="numeric",length=0)
for(i in 1:nrow(m1)){
m3<-m1[i,1]*m2
for(j in 2:ncol(m1))m3<-cbind(m3,m1[i,j]*m2)
if(i==1)kron<-m3
if(i>=2)kron<-rbind(kron,m3)
}
kron
}
