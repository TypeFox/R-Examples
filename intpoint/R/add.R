add <-
function(A,s){
if(is.array(A)){
r1<-nrow(A)
for(i in 1:r1){
aux<-array(0,c(r1,1))
aux[i]<-s
A<-cbind(A,aux)
}
}
else if(is.vector(A)){
l<-length(A)
aux<-array(0,c(1,l+1))
for(i in 1:l)
aux[i]<-A[i]
aux[l+1]<-s
A<-aux
}
return(A)
}
