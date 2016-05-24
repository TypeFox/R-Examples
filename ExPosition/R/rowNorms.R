rowNorms <- function(X,type=NULL,center=FALSE,scale=FALSE){
	
	if(is.null(type)){
		return(X)
	}else if(type=='hellinger'){
		#return((X/repmat(rowSums(X),1,ncol(X)))^(1/2))
		return( sqrt(X/repmat(rowSums(X),1,ncol(X))) )
	}else if(type == 'ca'){
		return(X/matrix(rowSums(X),nrow(X),ncol(X)))
	}else if (type == 'z'){
		return(t(expo.scale(t(X),center=TRUE,scale=TRUE)))
	}else if(type == 'other'){
		return(t(expo.scale(t(X),center=center,scale=scale)))
	}else{
		return(X)
	}
	
}