Cc <-
function(A){

A=as.matrix(A)
Ac=A-matrix(1,nrow(A),1)%*%SUM(A)$mc
return(Ac)
}
