truelambda <-
function(A,vector="n",accuracy=1e-7,iterations=1e+5){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
M<-A
eigvals<-eigen(M)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
A<-M/lambda
acr<-accuracy
if(!all(acr>=0&acr<1)) stop("accuracy must be between 0 and 1")
if(vector[1]=="n"){
    I<-diag(order)
    options(warn=-1)
    Nt<-project(A,time=10*order)
    options(warn=0)
    lambdat<-numeric(10*order)
    lambdat<-rep(list(lambdat),order)
    lastntminus1<-as.list(numeric(order))
    lastnt<-as.list(numeric(order))
    t<-numeric(order)
    truelamb<-matrix(0,ncol=2,nrow=order)
    for(i in 1:order){
        for(j in 1:((10*order)-1)){
            lambdat[[i]][j]<-Nt[j+1,i]/Nt[j,i]
        }
        lastntminus1[[i]]<-(A%^%((10*order)-1))%*%I[,i]
        lastnt[[i]]<-A%*%lastntminus1[[i]]
        if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
        t[i]<-1
        while(!all(lambdat[[i]]>=(lambdat[[i]][1]*(1-acr))&lambdat[[i]]<=(lambdat[[i]][1]*(1+acr)))&t[i]<iterations){
            lambdat[[i]][1:((10*order)-1)]<-lambdat[[i]][2:(10*order)]
            lastntminus1[[i]]<-A%*%lastntminus1[[i]]
            lastnt[[i]]<-A%*%lastnt[[i]]
            lambdat[[i]][10*order]<-sum(lastnt[[i]])/sum(lastntminus1[[i]])
            t[i]<-t[i]+1
            if(sum(lastnt[[i]])<.Machine$double.xmin|sum(lastnt[[i]])>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
        }
        if(!t[[i]]<iterations){
             stop("Model is not converging, try decreasing accuracy or increasing iterations\n(some imprimitive matrices may also not converge)")
        }
        truelamb[i,]<-c(lambdat[[i]][1]*lambda*(1-acr),lambdat[[i]][1]*lambda*(1+acr))
    }
}
else{
    n0<-vector
    vector<-n0/sum(n0)
    options(warn=-1)
    Nt<-project(A,vector=vector,standard.vec=TRUE,time=10*order)
    options(warn=0)
    lambdat<-numeric(10*order)
    for(i in 1:(10*order)){
        lambdat[i]<-Nt[i+1]/Nt[i]
    }
    lastntminus1<-(A%^%((10*order)-1))%*%vector
    lastnt<-A%*%lastntminus1
    if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
    t<-1
    while(!all(lambdat>=(lambdat[1]*(1-acr))&lambdat<=(lambdat[1]*(1+acr)))&t<iterations){
        lambdat[1:((10*order)-1)]<-lambdat[2:(10*order)]
        lastntminus1<-A%*%lastntminus1
        lastnt<-A%*%lastnt
        lambdat[10*order]<-sum(lastnt)/sum(lastntminus1)
        t<-t+1
        if(sum(lastnt)<.Machine$double.xmin|sum(lastnt)>.Machine$double.xmax) stop("projection calculation exceeded max or min normalized floating-point")
    }
    if(!t<iterations){
         stop("Model is not converging, try decreasing accuracy or increasing iterations\n(some imprimitive matrices may also not converge)")
    }
    truelamb<-c(lambdat[1]*lambda*(1-acr),lambdat[1]*lambda*(1+acr))
}
return(print(truelamb,-log10(acr)))
}

