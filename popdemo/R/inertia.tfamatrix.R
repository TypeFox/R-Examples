inertia.tfamatrix<-function(A, bound=NULL, vector="n", elementtype=NULL, Flim=c(-1,10), Plim=c(-1,10), plength=100, digits=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
laymat<-numeric(length(A))
dim(laymat)<-dim(A)
laymat[which(A>0)]<-1
laymat<-t(t(laymat)*cumsum(t(laymat)))
if(!is.null(elementtype)){
    type<-elementtype
}
if(is.null(elementtype)){
    type<-A
    type[type==0]<-NA
    type[1,][!is.na(type[1,])]<-"F"
    type[2:order,][!is.na(type[2:order,])&type[2:order,]<=1]<-"P"
    type[2:order,][!is.na(type[2:order,])&!type[2:order,]=="P"]<-"F"
}
p<-numeric(order*order*plength)
dim(p)<-c(order,order,plength)
lambda<-p
inertia<-p
for(i in 1:order){
    for(j in 1:order){
        if(A[i,j]!=0){
            d<-matrix(0,order)
            d[i,1]<-1
            e<-matrix(0,order)
            e[j,1]<-1
            if(type[i,j]=="P"){
                minp<-Plim[1]*A[i,j]
                maxp<-Plim[2]*A[i,j]-A[i,j]
                if(maxp>(1-A[i,j])) maxp<-(1-A[i,j])
                if(minp<(-A[i,j])) minp<-(-A[i,j])
                pert<-seq(minp,maxp,(maxp-minp)/(plength-1))
            }
            if(type[i,j]=="F"){
                minp<-Flim[1]*A[i,j]
                maxp<-Flim[2]*A[i,j]-A[i,j]
                if(minp<(-A[i,j])) minp<-(-A[i,j])
                pert<-seq(minp,maxp,(maxp-minp)/(plength-1))
            }
            transfer<-inertia.tfa(A,d,e,bound=bound,vector=vector,prange=pert,digits=digits)
            p[i,j,]<-c(transfer$p,rep(NA,plength-length(transfer$p)))
            lambda[i,j,]<-c(transfer$lambda,rep(NA,plength-length(transfer$lambda)))
            inertia[i,j,]<-c(transfer$inertia,rep(NA,plength-length(transfer$inertia)))
        }
    }
}
final<-list(p=p,lambda=lambda,inertia=inertia,layout=laymat)
class(final)<-c("tfamatrix","list")
return(final)
}
