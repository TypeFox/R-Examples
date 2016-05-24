inertia.tfa <-
function(A,d,e,vector="n",bound=NULL,prange,digits=1e-10){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) stop("Matrix A is imprimitive")
eigvals<-eigen(A)$values
lmax<-which.max(Re(eigvals))
lambda<-Re(eigvals[lmax])
if(any(is.null(d),is.null(e))) stop("please specify a perturbation structure using d and e")
d<-as.matrix(d)
e<-as.matrix(e)
mineigvals<-eigen(A+(min(prange)*d%*%t(e)))$values
maxeigvals<-eigen(A+(max(prange)*d%*%t(e)))$values
minlmax<-which.max(Re(mineigvals))
minlambda<-Re(mineigvals[minlmax])
maxlmax<-which.max(Re(maxeigvals))
maxlambda<-Re(maxeigvals[maxlmax])
lambdarange<-seq(minlambda,maxlambda,(maxlambda-minlambda)/(length(prange)-1))
lambdarange<-lambdarange[round(lambdarange,-log10(digits))!=round(lambda,-log10(digits))]
pvals<-1/tf(A,lambdarange,d,e)
I<-diag(order)
c<-as.matrix(rep(1,order))
bigd<-matrix(c(0,1),ncol=1)%x%d
bige<-matrix(c(1,0),ncol=1)%x%e
bigA2<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%I)
if(vector[1]=="n"){
    if(!any(bound=="upper",bound=="lower")) stop('Please specify bound="upper", bound="lower" or specify vector')
    vrange<-as.list(numeric(length(lambdarange)))
    mvrange<-numeric(length(lambdarange))
    vectorrange<-rep(list(as.matrix(numeric(order))),length(lambdarange))
    bigA1range<-as.list(numeric(length(lambdarange)))
    inertia<-numeric(length(lambdarange))
    if(bound=="upper"){
        for(i in 1:length(lambdarange)){
            vrange[[i]]<-abs(t(e)%*%solve((lambdarange[i]*I)-A))
            mvrange[i]<-which.max(vrange[[i]])
            vectorrange[[i]][mvrange[i]]<-1
            bigA1range[[i]]<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vectorrange[[i]]%*%t(c)))
            inertia[i]<-tf(bigA1range[[i]],lambdarange[i],bigd,bige)/tf(bigA2,lambdarange[i],bigd,bige)
        }
    }
    if(bound=="lower"){
        for(i in 1:length(lambdarange)){
            vrange[[i]]<-abs(t(e)%*%solve((lambdarange[i]*I)-A))
            mvrange[i]<-which.min(vrange[[i]])
            vectorrange[[i]][mvrange[i]]<-1
            bigA1range[[i]]<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vectorrange[[i]]%*%t(c)))
            inertia[i]<-tf(bigA1range[[i]],lambdarange[i],bigd,bige)/tf(bigA2,lambdarange[i],bigd,bige)
        }
    }
}
else{
    vector<-vector/sum(vector)
    bigA1<-(matrix(c(1,0,0,1),ncol=2,byrow=TRUE)%x%A)+(matrix(c(0,1,0,0),ncol=2,byrow=TRUE)%x%(vector%*%t(c)))
    inertia<-numeric(length(lambdarange))
    for(i in 1:length(lambdarange)){
        inertia[i]<-tf(bigA1,lambdarange[i],bigd,bige)/tf(bigA2,lambdarange[i],bigd,bige)
    }
}
final<-list(p=pvals,lambda=lambdarange,inertia=inertia)
class(final)<-c("tfa","list")
return(final)
}


