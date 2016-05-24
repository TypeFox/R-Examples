tfsens <-
function(A,d=NULL,e=NULL,startval=0.001,tolerance=1e-10,return.fit=FALSE,plot.fit=FALSE){
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
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
startlambda<-lambda+startval
f4.1<-tf(A,b=d,c=e,z=startlambda)^2
f5.1<-tf(A,b=d,c=e,z=startlambda,exp=-2)
S1<-f4.1/f5.1
eps<-startval/2
limlambda<-lambda+eps
f4.2<-tf(A,b=d,c=e,z=limlambda)^2
f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
S2<-f4.2/f5.2
lambdas<-c(startlambda,limlambda)
sensitivities<-c(S1,S2)
while(abs(S1-S2)>tolerance){
    f4.1<-f4.2
    f5.1<-f5.2
    S1<-S2
    eps<-eps/2
    limlambda<-lambda+eps
    f4.2<-tf(A,b=d,c=e,z=limlambda)^2
    f5.2<-tf(A,b=d,c=e,z=limlambda,exp=-2)
    S2<-f4.2/f5.2
    lambdas<-c(lambdas,limlambda)
    sensitivities<-c(sensitivities,S2)
}
if(plot.fit) plot(lambdas,sensitivities,xlab="lambda",ylab="sensitivity")
if(!return.fit) return(S2) else(return(list(sens=S2,lambda.fit=lambdas,sens.fit=sensitivities)))
}
