signA <-
function(A,i){
for(j in 1:ncol(A))
A[i,j]<-(-1)*A[i,j]
return(A)}
