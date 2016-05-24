a_artificial <-
function(A,b){
rows<-nrow(A)
cols<-ncol(A)
a<-matrix(0,rows,cols+1)
for(i in 1:rows)
for(j in 1:cols)
a[i,j]<-A[i,j]
for(i in 1:rows){
a[i,cols+1]<-b[i]
for(j in 1:cols)
a[i,cols+1]<-a[i,cols+1]-A[i,j]}
return(a)}
