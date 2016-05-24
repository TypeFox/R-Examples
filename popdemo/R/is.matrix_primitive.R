is.matrix_primitive <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
powermatrix<-A%^%((order^2)-(2*order)+2)
minval=min(powermatrix)
if(minval>0){
    return(TRUE)
}
else{
    return(FALSE)
}}
