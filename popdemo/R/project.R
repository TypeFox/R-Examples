project <-
function(A,vector="n",time=100,standard.A=FALSE,standard.vec=FALSE,return.vec=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)){
    warning("Matrix is reducible")
}
else{
    if(!is.matrix_primitive(A)) warning("Matrix is imprimitive")
}
if(standard.A==TRUE){
    M<-A
    eigvals<-eigen(M)$values
    lmax<-which.max(Re(eigvals))
    lambda<-Re(eigvals[lmax])
    A<-M/lambda
}
if(vector[1]=="n"){
    I<-diag(order)
    EmptyVec<-matrix(0,nrow=order,ncol=time+1)
    VecBias<-numeric((time+1)*order*order)
    dim(VecBias)<-c(time+1,order,order)
    dimnames(VecBias)[2]<-list(paste(rep("Stage",order),1:order,sep=""))
    dimnames(VecBias)[3]<-list(paste(rep("Bias",order),1:order,sep=""))
    PopBias<-numeric((time+1)*order)
    dim(PopBias)<-c(time+1,order)
    dimnames(PopBias)[2]<-list(paste(rep("Bias",order),1:order,sep=""))
    for (i in 1:order){
        VecBias[1,,i]<-I[,i]
        PopBias[1,i]<-1
            for (j in 1:time){
            VecBias[j+1,,i]<-A%*%VecBias[j,,i]
            PopBias[j+1,i]<-sum(VecBias[j+1,,i])
        }
    }
    class(PopBias)<-c("projection","matrix")
    if(return.vec){
        final<-list(N=PopBias,vec=VecBias)
        class(final)<-c("projection","list")
        return(final)
    }
    else{
        return(PopBias)
    }
}
else{
    n0<-vector
    if(standard.vec){
        vector<-n0/sum(n0)
    }
    Vec<-matrix(0,ncol=order,nrow=time+1)
    Vec[1,]<-vector
    dimnames(Vec)[2]<-list(paste(rep("Stage",order),1:order,sep=""))
    Pop<-numeric(time+1)
    Pop[1]<-sum(vector)
    for(i in 1:time){
        Vec[(i+1),]<-A%*%Vec[i,]
        Pop[i+1]<-sum(Vec[(i+1),])
    }
    class(Pop)<-c("projection","numeric")
    if(return.vec){
        final<-list(N=Pop,vec=Vec)
        class(final)<-c("projection","list")
        return(final)
    }
    else{
        return(Pop)
    }
}}

