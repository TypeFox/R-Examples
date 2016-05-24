dr <-
function(A,return.time=FALSE,x=10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
lambda2<-sort(Mod(eigvals),decreasing=TRUE)[2]
dr<-lambda/lambda2
if(return.time){
    t<-log(x)/log(dr)
    return(list(dr=dr,t=t))
}
else{
    return(dr)
}}
