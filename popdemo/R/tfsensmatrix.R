tfsensmatrix <- 
function(A,startval=0.001,tolerance=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("matrix is reducible")
if(!is.matrix_primitive(A)) warning("matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-eigvals[lmax]
if(max(Im(lambda))>0) stop("dominant eigenvalue contains nonzero imaginary component")
lambda<-Re(lambda)
if(tolerance>startval) stop("tolerance must be smaller than startval")
S<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            startlambda<-lambda+startval
            f4.1<-tf(A,b=d,c=e,z=startlambda)^2
            f5.1<-tf(A,b=d,c=e,z=startlambda,exp=-2)
            S1<-f4.1/f5.1
            eps<-startval/2
            limlambda<-lambda+eps
            f4.2<-tf(A,b=d,c=e,z=limlambda)^2
            f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
            S2<-f4.2/f5.2
            while(abs(S1-S2)>tolerance){
                f4.1<-f4.2
                f5.1<-f5.2
                S1<-S2
                eps<-eps/2
                limlambda<-lambda+eps
                f4.2<-tf(A,b=d,c=e,z=limlambda)^2
                f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
                S2<-f4.2/f5.2
            }
        S[i,j]<-S2
        }
    }
}
return(S)
}
