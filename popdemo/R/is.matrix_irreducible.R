is.matrix_irreducible <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
I<-diag(order)
IplusA<-I+A
powermatrix<-IplusA%^%(order-1)
minval<-min(powermatrix)
if(minval>0){
    return(TRUE)
}
else{
    return(FALSE)
}}

