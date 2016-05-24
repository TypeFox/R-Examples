tfa <-
function(A,d,e,prange=NULL,lambdarange=NULL,digits=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
if(is.null(lambdarange)&is.null(prange)) stop("Please specify one of either lambdarange or prange")
if(!is.null(lambdarange)&!is.null(prange)) stop("Please specify only one of either lambdarange or prange")
if(!is.null(lambdarange)&is.null(prange)){
    lambdarange<-lambdarange[round(lambdarange,-log10(digits))!=round(lambda,-log10(digits))]
    pvals<-1/tf(A,lambdarange,d,e)
}
if(is.null(lambdarange)&!is.null(prange)){
    mineigvals<-eigen(A+(min(prange)*d%*%t(e)))$values
    maxeigvals<-eigen(A+(max(prange)*d%*%t(e)))$values
    minlmax<-which.max(Re(mineigvals))
    minlambda<-Re(mineigvals[minlmax])
    maxlmax<-which.max(Re(maxeigvals))
    maxlambda<-Re(maxeigvals[maxlmax])
    lambdarange<-seq(minlambda,maxlambda,(maxlambda-minlambda)/(length(prange)-1))
    lambdarange<-lambdarange[round(lambdarange,-log10(digits))!=round(lambda,-log10(digits))]
    pvals<-1/tf(A,lambdarange,d,e)
}
final<-list(p=pvals,lambda=lambdarange)
class(final)<-c("tfa","list")
return(final)
}

