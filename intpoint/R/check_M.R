check_M <-
function(M,A){
if(is.array(M)){
r1<-nrow(A)+1
ca<-ncol(A)
cm<-ncol(M)
col<-ca-cm
rm<-nrow(M)
if(col!=0){
aux<-array(0,c(rm,col))
M<-cbind(M,aux)
}
A<-rbind(A,M)
r2<-nrow(A)
for(i in r1:r2){
aux<-array(0,c(r2,1))
aux[i]<--1
A<-cbind(A,aux)
}
}
else if(is.vector(M)){
r1<-nrow(A)+1
ca<-ncol(A)
cm<-length(M)
col<-ca-cm
aux1<-array(0,c(1,length(M)))
for(i in 1:length(M))
aux1[i]<-M[i]
M<-aux1
if(col!=0){
aux<-array(0,c(1,col))
M<-cbind(M,aux)
}
A<-rbind(A,M)
r2<-nrow(A)
for(i in r1:r2){
aux<-array(0,c(r2,1))
aux[i]<--1
A<-cbind(A,aux)
}

}
return(A)}
