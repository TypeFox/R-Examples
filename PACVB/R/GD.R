
.onLoad <- function(libname,pkgname){
	loadRcppModules()
}

GDHinge<-
function(X,Y,lambda,theta=0,K=100,v=10,ls=FALSE,B=0,family="F1",eps=0.05){
	skip=FALSE
	if(!is.matrix(X) | prod(dim(X)<1)){ warning("the design should be a matrix");skip=TRUE}
	if(is.vector(Y)) Y<-as.matrix(Y);
	if(ncol(Y)>1 & nrow(Y)==1) Y=t(Y)
	if(ncol(Y)!=1 & nrow(Y)!=1){ warning("response vector should be of size (n,1)");skip=TRUE}
	if(!prod(Y==-1| Y==1)){ warning("response vector should take values  {-1,1}");skip=TRUE}
	p<-ncol(X)
    if(length(theta)==1) theta<-rnorm(p+1)*0.01
    if(family=="F1"){ 
        fam=1
    }else{
        fam=0
    }
    if(skip)
	{
		return(1)
	}else{	
		return(GDHingeCxx(as.matrix(theta),X, Y,K,lambda,v,ls,B,fam,eps))
	}
}

GDHingeAUC<-
function (X,Y,lambda,theta=0,K=100,v=10,ls=FALSE,B=0,family="F1",eps=0.05){
	skip=FALSE
	if(!is.matrix(X) | prod(dim(X)<1)){ warning("the design should be a matrix");skip=TRUE}
	if(is.vector(Y)) Y<-as.matrix(Y);
	if(ncol(Y)>1 & nrow(Y)==1) Y=t(Y)
	if(ncol(Y)!=1 & nrow(Y)!=1){ warning("response vector should be of size (n,1)");skip=TRUE}
	if(!prod(Y==-1| Y==1)){ warning("response vector should take values  {-1,1}");skip=TRUE}
	p<-ncol(X)
    if(length(theta)==1) theta<-rnorm(p+1)*0.01
	if(family=="F1"){ fam=1
    }else{
        fam=0
    }
    if(skip)
	{
		return(1)
	}else{	
		return(GDHingeAUCCxx(as.matrix(theta),X, Y,K,lambda,v,ls,B,fam,eps))
	}
}
