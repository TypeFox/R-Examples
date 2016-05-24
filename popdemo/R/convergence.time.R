convergence.time <-
function(A,vector="n",accuracy=1e-2,iterations=1e+5){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix is reducible")
if(!is.matrix_primitive(A)) stop("Matrix is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
acr<-accuracy
if(!all(acr>=0&acr<1)) stop("accuracy must be between 0 and 1")
if(vector[1]=="n"){
    I<-diag(order)
    Nt<-project(A,time=10*order)
    lambdat<-numeric(10*order)
    lambdat<-rep(list(lambdat),order)
    lastntminus1<-as.list(numeric(order))
    lastnt<-as.list(numeric(order))
    t<-numeric(order)
    for(i in 1:order){
        for(j in 1:(10*order)){
            lambdat[[i]][j]<-Nt[j+1,i]/Nt[j,i]
        }
        lastntminus1[[i]]<-(A%^%((10*order)-1))%*%I[,i]
        lastnt[[i]]<-A%*%lastntminus1[[i]]
        t[i]<-1
        while(!all(lambdat[[i]]>(lambda*(1-acr))&lambdat[[i]]<(lambda*(1+acr)))&t[i]<iterations){
            lambdat[[i]][1:((10*order)-1)]<-lambdat[[i]][2:(10*order)]
            lastntminus1[[i]]<-A%*%lastntminus1[[i]]
            lastnt[[i]]<-A%*%lastnt[[i]]
            lambdat[[i]][10*order]<-sum(lastnt[[i]])/sum(lastntminus1[[i]])
            t[i]<-t[i]+1
            if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("Projection calculation exceeded max or min normalized floating-point")
        }
        if(!t[i]<iterations) stop("Model is not converging.  Reduce accuracy or increase iterations")
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    Nt<-project(A,vector=vector,time=10*order)
    lambdat<-numeric(10*order)
    for(i in 1:(10*order)){
        lambdat[i]<-Nt[i+1]/Nt[i]
    }
    lastntminus1<-(A%^%((10*order)-1))%*%vector
    lastnt<-A%*%lastntminus1
    t<-1
    while(!all(lambdat>(lambda*(1-acr))&lambdat<(lambda*(1+acr)))&t<iterations){
        lambdat[1:((10*order)-1)]<-lambdat[2:(10*order)]
        lastntminus1<-A%*%lastntminus1
        lastnt<-A%*%lastnt
        lambdat[10*order]<-sum(lastnt)/sum(lastntminus1)
        t<-t+1
        if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("Projection calculation exceeded max or min normalized floating-point")
    }
    if(!t<iterations) stop("Model is not converging.  Reduce accuracy or increase iterations")
}
return(t)}

