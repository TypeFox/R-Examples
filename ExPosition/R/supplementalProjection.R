supplementalProjection <- function(sup.transform=NULL,f.scores=NULL,Dv=NULL,scale.factor=NULL,symmetric=TRUE){
	if(is.null(sup.transform) || is.null(f.scores) || is.null(Dv)){
		stop('No inputs can be NULL.')
	}
	if(ncol(sup.transform)!=nrow(f.scores)){
		stop('Column dim of sup.transform does not match row dim of f.scores')
	}
	if(ncol(f.scores)!=length(Dv)){
		stop('Column dim of f.scores does not match length of Dv')
	}

	if(!symmetric){
		f.out <- sup.transform %*% f.scores	
	}else{
		f.out <- sup.transform %*% f.scores * matrix(Dv^-1,nrow(sup.transform),ncol(f.scores),byrow=TRUE)
	}
	if(!is.null(scale.factor)){		
		f.out <- f.out * matrix(scale.factor,nrow(f.out),ncol(f.scores),byrow=TRUE)
	}
	f.out <- replace(f.out,is.nan(f.out),0)
	d.out <- rowSums(f.out^2)
	r.out <- repmat((1/d.out),1,length(Dv)) * (f.out^2)
	r.out <- replace(r.out,is.nan(r.out),0)
	d.out <- as.matrix(d.out)
	return(list(f.out=f.out,d.out=d.out,r.out=r.out))
}