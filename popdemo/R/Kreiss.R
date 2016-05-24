Kreiss <-
function(A,bound=NULL,return.r=FALSE,theta=1,rlimit=100,step1=1e-3,step2=1e-6){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
I<-diag(order)
if(is.null(bound)) stop('Please specify either bound="upper or bound="lower"')
if(bound=="upper"){
    r1<-theta+step1
    K1<-(r1-theta)*(norm(solve((r1*I)-A)))
    Kreissbound1<-K1
    while(K1>=Kreissbound1 & r1<rlimit){
        Kreissbound1<-K1
        r1<-r1+step1
        K1<-(r1-theta)*(norm(solve((r1*I)-A)))
    }
    if(!r1<rlimit) stop("Maximum not reached and rlimit exceeded")
    r2<-r1
    K2<-(r2-theta)*(norm(solve((r2*I)-A)))
    Kreissbound2<-K1
    while(K2>=Kreissbound2 & r2>=(theta+(2*step2))){
        Kreissbound2<-K2
        r2<-r2-step2
        K2<-(r2-theta)*(norm(solve((r2*I)-A)))
    }
    rstar<-r2+step2
    Kstar<-Kreissbound2
    r.1<-r2
    K.1<-K2
    if(r2<(theta+(2*step2))){
        Kreissbound<-K.1
        r<-r.1
    }
    else{
        Kreissbound<-Kstar
        r<-rstar
    }
    if(return.r) return(list(Kreissbound=Kreissbound,r=r)) else(return(Kreissbound))
}
if(bound=="lower"){
    r1<-theta+step1
    K1<-(r1-theta)*(minCS(solve((r1*I)-A)))
    Kreissbound1<-K1
    while(K1<=Kreissbound1 & r1<rlimit){
        Kreissbound1<-K1
        r1<-r1+step1
        K1<-(r1-theta)*(minCS(solve((r1*I)-A)))
    }
    if(!r1<rlimit) stop("Minimum not reached and rlimit exceeded")
    r2<-r1
    K2<-(r2-theta)*(minCS(solve((r2*I)-A)))
    Kreissbound2<-K1
    while(K2<=Kreissbound2 & r2>=(theta+(2*step2))){
        Kreissbound2<-K2
        r2<-r2-step2
        K2<-(r2-theta)*(minCS(solve((r2*I)-A)))
    }
    rstar<-r2+step2
    Kstar<-Kreissbound2
    r.1<-r2
    K.1<-K2
    if(r2<(theta+(2*step2))){
        Kreissbound<-K.1
        r<-r.1
    }
    else{
        Kreissbound<-Kstar
        r<-rstar
    }
    if(return.r) return(list(Kreissbound=Kreissbound,r=r)) else (return(Kreissbound))
}}

