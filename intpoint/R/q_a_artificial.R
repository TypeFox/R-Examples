q_a_artificial <-
function(A,j){
rows<-nrow(A)
for(i in 1:rows)
A[i,j]<-0
return(A)}
