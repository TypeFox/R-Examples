check_mm <-
function(m,A){
if(is.vector(m)||is.array(m)){
r1<-nrow(A)+1
A<-rbind(A,m)
r2<-nrow(A)
for(i in r1:r2){
aux<-array(0,c(r2,1))
aux[i]<-1
A<-cbind(A,aux)
}
}
return(A)}
