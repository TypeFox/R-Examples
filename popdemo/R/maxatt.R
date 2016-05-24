maxatt <-
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
    minN<-numeric(order)
    times<-numeric(order)
    for(i in 1:order){
        minN[i]<-min(projection[,i])
        times[i]<-which.min(projection[,i])-1
    }
    rhomin<-min(minN)
    t<-times[which.min(minN)]
    stage<-which.min(minN)
    if(return.N){
        Nt<-rhomin*lambda^t
        if(all(return.t,return.stage)) return(list(maxatt=rhomin,N=Nt,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxatt=rhomin,N=Nt,t=t))
        if(all(!return.t,return.stage)) return(list(maxatt=rhomin,N=Nt,stage=stage))
        if(!all(return.t,return.stage)) return(list(maxatt=rhomin,N=Nt))
    }
    else{
        if(all(return.t,return.stage)) return(list(maxatt=rhomin,t=t,stage=stage))
        if(all(return.t,!return.stage)) return(list(maxatt=rhomin,t=t))
        if(all(!return.t,return.stage)) return(list(maxatt=rhomin,stage=stage))
        if(!all(return.t,return.stage)) return(rhomin)
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    maxtime<-convergence.time(A,vector=vector,accuracy=conv.accuracy,iterations=conv.iterations)
    projection<-project(A,vector=vector,time=maxtime)
    rhomin<-min(projection)
    t<-which.min(projection)-1
    if(rhomin<1){
        if(return.N){
            Nt<-rhomin*sum(n0)*lambda^t
            if(return.t) return(list(maxatt=rhomin,N=Nt,t=t)) else(return(list(maxatt=rhomin,N=Nt)))
        }
        else{
            if(return.t) return(list(maxatt=rhomin,t=t)) else(return(rhomin))
        }
    }
    else{
        stop("Model does not attenuate.  Cannot compute maximum attenuation")
    }
}}

