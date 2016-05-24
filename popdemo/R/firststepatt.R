firststepatt <-
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
    firstep<-minCS(A)
    if(return.N){
        Nt<-firstep*lambda
        return(list(firststepatt=firstep,N=Nt))
    }
    else{
        return(firstep)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    firstep<-sum(A%*%vector)
    if(firstep<1){
        if(return.N){
            Nt<-firstep*sum(n0)*lambda
            return(list(firststepatt=firstep,N=Nt))
        }
        else{
            return(firstep)
        }
    }
    else{
        stop("Model amplifies.  Cannot compute first-timestep attenuation,\n use function reactivity to calculate reactivity")
    }
}}

