inertia <-
function(A,vector="n",bound=NULL,return.N=FALSE,t=NULL){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) stop("Matrix A is imprimitive")
M<-A
reigs<-eigen(A)
leigs<-eigen(t(A))
lmax<-which.max(Re(reigs$values))
lambda<-Re(leigs$values[lmax])
A<-M/lambda
w<-as.matrix(abs(Re(reigs$vectors[,lmax])))
v<-as.matrix(abs(Re(leigs$vectors[,lmax])))
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    if(bound=="upper"){
        rhoinfinity<-as.vector((max(v)*sum(w))/(t(v)%*%w))
        if(return.N){
            if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
            warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
            N<-rhoinfinity*lambda^t
            return(list(upper.inertia=rhoinfinity,N=N))
        }
        else{
            return(rhoinfinity)
        }
    }
    if(bound=="lower"){
        rhoinfinity<-as.vector((min(v)*sum(w))/(t(v)%*%w))
        if(return.N){
            if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
            warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
            N<-rhoinfinity*lambda^t
            return(list(lower.inertia=rhoinfinity,N=N))
        }
        else{
            return(rhoinfinity)
        }
    }
}
else{
    if(!is.null(bound)) warning("Specification of vector overrides calculation of bound")
    n0<-vector
    vector<-n0/sum(n0)
    Pinfinity<-as.vector((t(v)%*%vector*sum(w))/(t(v)%*%w))
    if(return.N){
        if(is.null(t)) stop("Please specify a value of t at which N is to be calculated")
        warning("Estimation of N will be  inaccurate for\n t where the model has not converged.")
        N<-Pinfinity*sum(n0)*lambda^t
        return(list(inertia=Pinfinity,N=N))
    }
    else{
        return(Pinfinity)
    }
}}

