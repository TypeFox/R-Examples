maxamp <-
function(A,vector="n",return.N=FALSE,return.t=FALSE,return.stage=FALSE,conv.iterations=1e+5,conv.accuracy=1e-5){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
if(vector[1]=="n"){
    maxtime<-max(convergence.time(A,accuracy=conv.accuracy,iterations=conv.iterations))
    projection<-project(A,time=maxtime)
    maxN<-numeric(order)
    times<-numeric(order)
    for(i in 1:order){
        maxN[i]<-max(projection[,i])
        times[i]<-which.max(projection[,i])-1
    }
    rhomax<-max(maxN)
    t<-times[which.max(maxN)]
    stage<-which.max(maxN)
    if(return.N){
        Nt<-rhomax*lambda^t
        if(all(return.t,return.stage)) return(list(maxamp=rhomax,N=Nt,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxamp=rhomax,N=Nt,t=t))
        if(all(!return.t,return.stage)) return(list(maxamp=rhomax,N=Nt,stage=stage))
        if(!all(return.t,return.stage)) return(list(maxamp=rhomax,N=Nt))
    }
    else{
        if(all(return.t,return.stage)) return(list(maxamp=rhomax,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxamp=rhomax,t=t))
        if(all(!return.t,return.stage)) return(list(maxamp=rhomax,stage=stage))
        if(!all(return.t,return.stage)) return(rhomax)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    maxtime<-convergence.time(A,vector=vector,accuracy=conv.accuracy,iterations=conv.iterations)
    projection<-project(A,vector=vector,time=maxtime)
    rhomax<-max(projection)
    t<-which.max(projection)-1
    if(rhomax>1){
        if(return.N){
            Nt<-rhomax*sum(n0)*lambda^t
            if(return.t) return(list(maxamp=rhomax,N=Nt,t=t)) else(return(list(maxamp=rhomax,N=Nt)))
        }
        else{
            if(return.t) return(list(maxamp=rhomax,t=t)) else(return(rhomax))
        }
    }
    else{
        stop("Model does not amplify.  Cannot compute maximum amplification")
    }
}}

