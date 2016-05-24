
.onLoad <- function(libname,pkgname){
	loadRcppModules()
}

EPprobit <-
function (X, Y, s){
	skip=FALSE
	if(s<=0){ warning("variance should be positive");skip=TRUE}
	if(!is.matrix(X) | prod(dim(X)<1)){ warning("the design should be a matrix");skip=TRUE}
	if(is.vector(Y)) Y<-as.matrix(Y);
	if(ncol(Y)>1 & nrow(Y)==1) Y=t(Y)
	if(ncol(Y)!=1 & nrow(Y)!=1){ warning("response vector should be of size (n,1)");skip=TRUE}
	if(!prod(Y==0 | Y==1)){ warning("response vector should take values  {0,1}");skip=TRUE}
	if(skip)
	{
		return(1)
	}else{	
		return(EPprobitCxx(X, Y, s))
	}
}

EPlogit <-
function (X, Y, s){ 
	skip=FALSE
	if(s<=0){ warning("variance should be positive");skip=TRUE}
	if(!is.matrix(X) | prod(dim(X)<1)){ warning("the design should be a matrix");skip=TRUE}
	if(is.vector(Y)) Y<-as.matrix(Y);
	if(ncol(Y)>1 & nrow(Y)==1) Y=t(Y)
	if(ncol(Y)!=1 & nrow(Y)!=1){ warning("response vector should be of size (n,1)");skip=TRUE}
	if(!prod(Y==0 | Y==1)){ warning("response vector should take values {0,1}");skip=TRUE}
	if(skip)
	{
		return(1)
	}else{	
		return(EPlogitCxx(X, Y, s))
	}
}
