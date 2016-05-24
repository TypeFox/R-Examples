reactivity <-
function(A,vector="n",return.N=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)){
    warning("Matrix is reducible")
}
else{
    if(!is.matrix_primitive(A)) warning("Matrix is imprimitive")
}
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
if(vector[1]=="n"){
    reac<-norm(A)
    if(return.N){
        N<-reac*lambda
        return(list(reactivity=reac,N=N))
    }
    else{
        return(reac)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    reac<-sum(A%*%vector)
    if(reac>1){
        if(return.N){
            N<-reac*sum(n0)*lambda
            return(list(reactivity=reac,N=N))
        }
        else{
            return(reac)
        }
    }
    else{
        stop("Model attenuates.  Cannot compute reactivity, use function\n firststepatt to calculate first-timestep attenuation")
    }
}}

