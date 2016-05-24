inertia.tfsensmatrix <-
function(A,vector="n",bound=NULL,startval=0.001,tolerance=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix is reducible")
if(!is.matrix_primitive(A)) stop("Matrix is imprimitive")
leigs<-eigen(t(A))
lmax<-which.max(Re(leigs$values))
lambda<-Re(leigs$values[lmax])
v<-as.matrix(abs(Re(leigs$vectors[,lmax])))
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    n0<-as.matrix(rep(0,order))
    if(bound=="upper") n0[which.max(v),1]<-1
    if(bound=="lower") n0[which.min(v),1]<-1
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-as.matrix(vector/sum(vector))
}
if(tolerance>startval) stop("tolerance must be smaller than startval")
ones<-as.matrix(rep(1,order))
S<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            startlambda<-lambda+startval
            f1.1<-tf(A,b=n0,c=e,z=startlambda)
            f1d.1<--tf(A,b=n0,c=e,z=startlambda,exp=-2)
            f2.1<-tf(A,b=d,c=ones,z=startlambda)
            f2d.1<--tf(A,b=d,c=ones,z=startlambda,exp=-2)
            f3.1<-tf(A,b=d,c=e,z=startlambda,exp=-2)
            f3d.1<--2*tf(A,b=d,c=e,z=startlambda,exp=-3)
            f4.1<-tf(A,b=d,c=e,z=startlambda)^2
            f5.1<-tf(A,b=d,c=e,z=startlambda,exp=-2)
            diff1.1<-((f1d.1*f2.1*f3.1)+(f1.1*f2d.1*f3.1)-(f1.1*f2.1*f3d.1))/f3.1^2
            diff2.1<-f4.1/f5.1
            S1<-diff1.1*diff2.1
            eps<-startval/2
            limlambda<-lambda+eps
            f1.2<-tf(A,b=n0,c=e,z=limlambda)
            f1d.2<--tf(A,b=n0,c=e,z=limlambda,exp=-2)
            f2.2<-tf(A,b=d,c=ones,z=limlambda)
            f2d.2<--tf(A,b=d,c=ones,z=limlambda,exp=-2)
            f3.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
            f3d.2<--2*tf(A,b=d,c=e,z=limlambda,exp=-3)
            f4.2<-tf(A,b=d,c=e,z=limlambda)^2
            f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
            diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
            diff2.2<-f4.2/f5.2
            S2<-diff1.2*diff2.2
            while(abs(S1-S2)>tolerance){
                f1.1<-f1.2
                f1d.1<-f1d.2
                f2.1<-f2.2
                f2d.1<-f2d.2
                f3.1<-f3.2
                f3d.1<-f3d.2
                f4.1<-f4.2
                f5.1<-f5.2
                diff1.1<-diff1.2
                diff2.1<-diff2.2
                S1<-S2
                eps<-eps/2
                limlambda<-lambda+eps
                f1.2<-tf(A,b=n0,c=e,z=limlambda)
                f1d.2<--tf(A,b=n0,c=e,z=limlambda,exp=-2)
                f2.2<-tf(A,b=d,c=ones,z=limlambda)
                f2d.2<--tf(A,b=d,c=ones,z=limlambda,exp=-2)
                f3.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
                f3d.2<--2*tf(A,b=d,c=e,z=limlambda,exp=-3)
                f4.2<-tf(A,b=d,c=e,z=limlambda)^2
                f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
                diff1.2<-((f1d.2*f2.2*f3.2)+(f1.2*f2d.2*f3.2)-(f1.2*f2.2*f3d.2))/f3.2^2
                diff2.2<-f4.2/f5.2
                S2<-diff1.2*diff2.2
            }
            S[i,j]<-S2
        }
    }
}
return(S)
}
